import gzip
from enum import Enum
from io import BufferedIOBase
from typing import Optional

from . import _bgzf, _zstd


class Compression(Enum):
    AUTO = "auto"
    NONE = "none"
    BGZF = "bgzf"
    GZIP = "gzip"
    ZSTD = "zstd"


def _format_enabled(compression: Compression) -> bool:
    if (compression is Compression.ZSTD) and (not _zstd._HAS_ZSTD):
        return False

    return True


def _is_format(source: BufferedIOBase, compression: Compression) -> bool:
    if compression is Compression.BGZF:
        return _bgzf.is_bgzf(source)
    elif compression is Compression.GZIP:
        return _bgzf.is_gzip(source)
    elif compression is Compression.ZSTD:
        return _zstd.is_zstd(source)


def _read_chunk(source: BufferedIOBase, compression: Compression) -> bytes:
    if compression is Compression.BGZF:
        return _bgzf._read_block(source)
    elif compression is Compression.ZSTD:
        return _zstd._read_frame(source)

    raise ValueError(f"Cannot read a chunk from compression type '{compression}'")


def _decompress(source: bytes, compression: Compression) -> bytes:
    if (compression is Compression.GZIP) or (compression is Compression.BGZF):
        return gzip.decompress(source)
    elif compression is Compression.ZSTD:
        return _zstd._decompress(source)

    raise ValueError(f"Cannot decompress from compression type '{compression}'")


def _index(
    source: BufferedIOBase, destination: BufferedIOBase, compression: Compression
) -> None:
    if compression is Compression.BGZF:
        return _bgzf.index(source, destination)
    elif compression is Compression.ZSTD:
        return _zstd.index(source, destination)

    raise ValueError(f"Could not create GZI index for compression '{compression}'")


def _detect_format(source: BufferedIOBase):
    if _bgzf.is_bgzf(source):
        return Compression.BGZF
    elif _bgzf.is_gzip(source):
        return Compression.GZIP
    elif _zstd.is_zstd(source):
        return Compression.ZSTD

    raise RuntimeError("Failed to detect compression type from source")


def _get_line_iterator(
    source: BufferedIOBase, compression: Optional[Compression] = None
):
    if compression is None:
        compression = _detect_format(source)

    if (compression is Compression.GZIP) or (compression is Compression.BGZF):
        return gzip.open(source, mode="rt")
    elif compression is Compression.ZSTD:
        return _zstd._get_line_iterator(source)

    raise ValueError(
        f"Cannot create line iterator from compression type '{compression}'"
    )
