import logging
import struct
import zlib
from collections import namedtuple
from io import BufferedIOBase
from typing import Optional

MAX_BLOCK_SIZE_IN_BYTES = 65_536
# solved for source_len using
# MAX_BLOCK_SIZE_IN_BYTES = _deflate_bound(
#   source_len, wrap_bytes = (20 + 6),
#   w_bits = 15, mem_level = 8)
UNCOMPRESSED_BLOCK_DATA_BOUND_IN_BYTES = 65_485
GZIP_HEADER = b"\x1f\x8b"
BGZF_HEADER = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
BGZF_HEADER_STRUCT_FORMAT = "<BBBBLBBHBBHH"
BGZF_TAILER_STRUCT_FORMAT = "<LL"


def _write_bgzf_block(destination: BufferedIOBase, content: bytes):
    if len(content) > UNCOMPRESSED_BLOCK_DATA_BOUND_IN_BYTES:
        raise ValueError(
            "Attempting to compress too much data into a single BGZF block"
        )

    # compress the source bytes
    compressor = zlib.compressobj(wbits=-15, memLevel=8)
    compressed = compressor.compress(content)
    compressed += compressor.flush(zlib.Z_FINISH)
    bsize = len(compressed) + 25

    if bsize >= MAX_BLOCK_SIZE_IN_BYTES:
        raise RuntimeError(
            "Data compressed too large to be fit into a single BGZF block"
        )

    # write ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN, SI1, SI2, SLEN
    destination.write(BGZF_HEADER)
    # write BSIZE
    destination.write(struct.pack("<H", bsize))
    # write CDATA
    destination.write(compressed)
    # write CRC32 and ISIZE
    destination.write(
        struct.pack(BGZF_TAILER_STRUCT_FORMAT, zlib.crc32(content), len(content))
    )


def compress(source: BufferedIOBase, destination: BufferedIOBase):
    block = source.read(UNCOMPRESSED_BLOCK_DATA_BOUND_IN_BYTES)
    while block:
        # try to split blocks on newlines if possible
        newline_idx = block.rfind(b"\n")

        if newline_idx >= 0:
            remainder = block[newline_idx + 1 :]
            block = block[: newline_idx + 1]
        _write_bgzf_block(destination, block)

        block = remainder + source.read(
            UNCOMPRESSED_BLOCK_DATA_BOUND_IN_BYTES - len(remainder)
        )
        remainder = b""

    # write EOF marker block
    _write_bgzf_block(destination, b"")


BGZFHeader = namedtuple(
    "BGZFHeader",
    [
        "ID1",
        "ID2",
        "CM",
        "FLG",
        "MTIME",
        "XFL",
        "OS",
        "XLEN",
        "SI1",
        "SI2",
        "SLEN",
        "BSIZE",
    ],
)
BGZFTailer = namedtuple("BGZFTailer", ["CRC32", "ISIZE"])


def _read_block(source: BufferedIOBase) -> bytes:
    block_header_length = struct.calcsize(BGZF_HEADER_STRUCT_FORMAT)
    block_tailer_length = struct.calcsize(BGZF_TAILER_STRUCT_FORMAT)

    block_header = BGZFHeader(
        *struct.unpack(BGZF_HEADER_STRUCT_FORMAT, source.read(block_header_length))
    )
    block_remainder = source.read(block_header.BSIZE + 1 - block_header_length)
    payload_compressed, block_tailer = block_remainder[
        :-block_tailer_length
    ], BGZFTailer(
        *struct.unpack(
            BGZF_TAILER_STRUCT_FORMAT, block_remainder[-block_tailer_length:]
        )
    )
    payload = zlib.decompress(payload_compressed, wbits=-zlib.MAX_WBITS)

    if zlib.crc32(payload) != block_tailer.CRC32:
        logging.warning("BGZF block CRC32 failed to validate")
    if len(payload) != block_tailer.ISIZE:
        logging.warning("BGZF block data size does not match metadata")

    return payload


def index(source: BufferedIOBase, destination: BufferedIOBase):
    block_header_length = struct.calcsize(BGZF_HEADER_STRUCT_FORMAT)
    block_tailer_length = struct.calcsize(BGZF_TAILER_STRUCT_FORMAT)

    blocks = []
    source.seek(0)
    compressed_offset = 0
    uncompressed_offset = 0
    block_header = source.read(block_header_length)
    while block_header:
        block_header = BGZFHeader(
            *struct.unpack(BGZF_HEADER_STRUCT_FORMAT, block_header)
        )
        blocks.append((compressed_offset, uncompressed_offset))
        # skip content
        source.seek(compressed_offset + block_header.BSIZE + 1 - block_tailer_length)
        block_tailer = BGZFTailer(
            *struct.unpack(BGZF_TAILER_STRUCT_FORMAT, source.read(block_tailer_length))
        )
        uncompressed_offset += block_tailer.ISIZE

        # set up next loop
        compressed_offset = source.tell()
        block_header = source.read(block_header_length)

    # don't include first and last (empty) block
    blocks = blocks[1:-1]
    destination.write(struct.pack("<Q", len(blocks)))
    for block in blocks:
        destination.write(struct.pack("<QQ", *block))


def is_gzip(stream: BufferedIOBase) -> bool:
    header = stream.read(len(GZIP_HEADER))
    # reset the offset
    stream.seek(stream.tell() - len(GZIP_HEADER))
    return header == GZIP_HEADER


def is_bgzf(stream: BufferedIOBase) -> bool:
    header = stream.read(len(BGZF_HEADER))
    # reset the offset
    stream.seek(stream.tell() - len(BGZF_HEADER))
    return header == BGZF_HEADER


def _deflate_bound(
    source_len: int,
    wrap_bytes: int = 0,
    w_bits: Optional[int] = None,
    mem_level: Optional[int] = None,
) -> int:
    """Computes an upper bound for the compressed size of a source with a known length.

    Implementation is modified from zlib.
    """

    fixed_len = (
        source_len + (source_len >> 3) + (source_len >> 8) + (source_len >> 9) + 4
    )
    store_len = (
        source_len + (source_len >> 5) + (source_len >> 7) + (source_len >> 11) + 7
    )

    if (w_bits is None) or (mem_level is None):
        return max(fixed_len, store_len) + wrap_bytes

    hash_bits = mem_level + 7
    if (w_bits != 15) or (hash_bits != 15):
        return (fixed_len if w_bits <= hash_bits else store_len) + wrap_bytes

    return (
        source_len
        + (source_len >> 12)
        + (source_len >> 14)
        + (source_len >> 25)
        + 13
        - 6
        + wrap_bytes
    )
