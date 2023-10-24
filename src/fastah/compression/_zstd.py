import struct
from io import BufferedIOBase, TextIOWrapper
from typing import Callable

try:
    import zstandard

    _HAS_ZSTD = True
except ModuleNotFoundError:
    _HAS_ZSTD = False


def _requires_zstandard(func: Callable) -> Callable:
    def wrapper(*args, **kwargs):
        if not _HAS_ZSTD:
            raise RuntimeError("zstandard must be installed")
        return func(*args, **kwargs)

    return wrapper


@_requires_zstandard
def is_zstd(stream: BufferedIOBase) -> bool:
    header = stream.read(4)
    # reset the offset
    stream.seek(stream.tell() - 4)
    return header == zstandard.FRAME_HEADER


@_requires_zstandard
def index(source: BufferedIOBase, destination: BufferedIOBase):
    blocks = []
    source.seek(0)
    compressed_offset = 0
    uncompressed_offset = 0
    compressor = zstandard.ZstdDecompressor()
    reader = compressor.stream_reader(source=source, read_across_frames=False)

    block = reader.read()
    while block:
        blocks.append((compressed_offset, uncompressed_offset))
        compressed_offset = source.tell()
        uncompressed_offset += len(block)
        block = reader.read()

    # don't include first block
    blocks = blocks[1:]
    destination.write(struct.pack("<Q", len(blocks)))
    for block in blocks:
        destination.write(struct.pack("<QQ", *block))


@_requires_zstandard
def _decompress(source: bytes) -> bytes:
    decompressor = zstandard.ZstdDecompressor()
    return decompressor.decompress(source)


@_requires_zstandard
def _read_frame(source: BufferedIOBase) -> bytes:
    decompressor = zstandard.ZstdDecompressor()
    reader = decompressor.stream_reader(source=source, read_across_frames=False)

    return reader.read()


@_requires_zstandard
def _get_line_iterator(source: BufferedIOBase) -> TextIOWrapper:
    decompressor = zstandard.ZstdDecompressor()
    byte_stream = decompressor.stream_reader(source=source)
    text_stream = TextIOWrapper(byte_stream, encoding="utf-8")

    return text_stream
