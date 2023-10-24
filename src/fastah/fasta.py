import gzip
import multiprocessing
import struct
from bisect import bisect_left
from collections import OrderedDict
from dataclasses import dataclass
from io import BufferedIOBase, IOBase, TextIOBase
from operator import attrgetter
from os import PathLike
from pathlib import Path
from typing import Generator, Optional, Union

from .compression import (
    Compression,
    _bgzf,
    _decompress,
    _format_enabled,
    _get_line_iterator,
    _is_format,
    _read_chunk,
)


class FASTARecord:
    _file: Optional["FASTAFile"] = None

    def __init__(self, id: str, sequence: str, description: str = ""):
        self._id = id
        self._sequence = sequence
        self._description = description

    @property
    def id(self) -> str:
        return self._id

    @property
    def description(self) -> str:
        return self._description

    @property
    def sequence(self) -> str:
        if self._file is not None:
            # cache fetched sequence
            self._sequence = self[:]
            # unlink from FASTAFile
            self._file = None

        return self._sequence

    def __len__(self) -> int:
        if self._file is not None:
            return self._file._get_seqid_length(self.id)
        else:
            return len(self.sequence)

    def __getitem__(self, key) -> str:
        if isinstance(key, slice):
            start = key.start
            stop = key.stop
            if key.step == 0:
                raise ValueError("slice step cannot be zero")
            step = key.step
        elif isinstance(key, int):
            if (key >= len(self)) or (key < -len(self)):
                raise IndexError("FASTA sequence index out of range")

            start = key
            stop = key + 1
            step = None
        else:
            raise ValueError(
                f"Only integers and slices can be used to index into a FASTA record, not '{type(key).__name__}'"
            )

        if self._sequence:
            return self._sequence[start:stop]
        else:
            return self._file._fetch(seqid=self.id, start=start, stop=stop, step=step)

    def __str__(self) -> str:
        return self.sequence

    def __repr__(self) -> str:
        if self._file is not None:
            sequence = "..."
        else:
            sequence = self.sequence if len(self) <= 3 else f"{self.sequence[:3]}..."

        return f"{self.__class__.__name__}(id='{self.id}', description='{self.description}', sequence='{sequence}')"

    def __getstate__(self):
        if self._file is not None:
            # cache sequence and unlink from backing FASTAFile
            _ = self.sequence

        return self.__dict__

    @staticmethod
    def _from_FASTA_file(file: "FASTAFile", seqid: str) -> "FASTARecord":
        record = FASTARecord(id=seqid, sequence="")
        record._file = file

        return record


@dataclass(frozen=True)
class FaiEntry:
    name: str
    length: int
    offset: int
    linebases: int
    linewidth: int
    qualoffset: Optional[int] = None


@dataclass(frozen=True)
class CompressedBlock:
    compressed_offset: int
    uncompressed_offset: int


def index(
    source: Union[PathLike, TextIOBase, BufferedIOBase, str],
    destination: Union[PathLike, TextIOBase, str],
):
    should_close_source = False
    if isinstance(source, PathLike) or isinstance(source, str):
        source = Path(source)
        source = open(source, mode="r" if (source.suffix == "fa") else "rb")
        should_close_source = True

    if isinstance(source, BufferedIOBase):
        try:
            source = _get_line_iterator(source)
        except RuntimeError as e:
            raise ValueError("Unknown compression format, cannot index") from e

    should_close_destination = False
    if isinstance(destination, PathLike) or isinstance(destination, str):
        destination = open(destination, mode="w")
        should_close_destination = True

    offset = 0
    seqid = None
    # Python does not allow getting the current file offset with
    # .tell() inside a for loop over the file, so use a while loop
    line = source.readline()
    encountered_blank_line_before = False
    while line:
        offset += len(line)

        line_rstripped = line.rstrip()

        if encountered_blank_line_before:
            raise RuntimeError(
                "Encountered a blank line in the middle of the FASTA file"
            )

        if not line_rstripped:
            encountered_blank_line_before = True
        elif line.startswith(">"):
            if seqid is not None:
                destination.write(
                    f"{seqid}\t{seq_len}\t{seq_offset}\t{seq_linebases}\t{seq_linewidth}\n"
                )
            seqid = line.rstrip().split(" ")[0][1:]
            seq_len = 0
            seq_offset = offset
            seq_linebases = None
            seq_linewidth = None
            encountered_unequal_linebases_before = False
            encountered_unequal_linewidths_before = False
        else:
            if seq_linebases is None:
                seq_linebases = len(line_rstripped)
            elif len(line_rstripped) != seq_linebases:
                if encountered_unequal_linebases_before:
                    raise RuntimeError(
                        f"Encountered unequal numbers of bases in lines of sequence record '{seqid}'"
                    )
                encountered_unequal_linebases_before = True

            if seq_linewidth is None:
                seq_linewidth = len(line)
            elif len(line) != seq_linewidth:
                terminator_width = len(line) - len(line_rstripped)
                if terminator_width != (seq_linewidth - seq_linebases):
                    # always error if unequal terminator widths
                    raise RuntimeError(
                        f"Encountered lines with unequal terminator widths in sequence record '{seqid}'"
                    )
                elif encountered_unequal_linewidths_before:
                    raise RuntimeError(
                        f"Encountered unequal line widths in sequence record '{seqid}'"
                    )
                encountered_unequal_linewidths_before = True

            seq_len += len(line_rstripped)

        line = source.readline()

    if seqid is not None:
        # write out final record index
        destination.write(
            f"{seqid}\t{seq_len}\t{seq_offset}\t{seq_linebases}\t{seq_linewidth}\n"
        )

    if should_close_source:
        source.close()

    if should_close_destination:
        destination.close()


class FASTAFile:
    _path: Optional[Path]
    _index: Optional[OrderedDict[str, FaiEntry]]
    _index_compressed: Optional[tuple[CompressedBlock, ...]]

    def __init__(
        self,
        source: Union[PathLike, TextIOBase, BufferedIOBase, str],
        index: Optional[Union[bool, PathLike, TextIOBase, str]] = True,
        index_compressed: Optional[Union[bool, PathLike, BufferedIOBase, str]] = True,
        compression: Union[Compression, str] = Compression.AUTO,
    ):
        # cast string arguments to appropriate types
        if isinstance(source, str):
            source = Path(source)
        if isinstance(index, str):
            index = Path(index)
        if isinstance(index_compressed, str):
            index_compressed = Path(index_compressed)
        if isinstance(compression, str):
            compression = Compression(compression)

        if isinstance(source, PathLike):
            source = Path(source)
            self._path = source
        else:
            self._path = None

        # detect compression
        if compression is Compression.AUTO:
            if isinstance(source, TextIOBase):
                # text streams won't be compressed
                compression = Compression.NONE
            elif isinstance(source, Path) and (source.suffix == ".fa"):
                compression = Compression.NONE

        self._stream = (
            source
            if isinstance(source, IOBase)
            else open(source, "r" if compression is Compression.NONE else "rb")
        )
        self._stream_lock = multiprocessing.Lock()

        # detect compression type if not uncompressed
        if compression is Compression.AUTO:
            # peek in binary stream for BGZF header
            if _is_format(self._stream, compression.BGZF):
                compression = Compression.BGZF
            elif _is_format(self._stream, compression.GZIP):
                compression = Compression.GZIP
                self._stream = gzip.open(self._stream, mode="rt")
            elif _format_enabled(compression.ZSTD) and _is_format(
                self._stream, compression.ZSTD
            ):
                compression = Compression.ZSTD
            else:
                raise ValueError(
                    "Failed to auto-detect FASTA compression type from provided inputs. If you are using zstd, make sure the zstandard package is installed."
                )

        self._compression = compression

        if index and isinstance(source, PathLike):
            # attempt to auto-locate the index from the source path
            index = source.with_suffix(source.suffix + ".fai")
        elif isinstance(index, PathLike) and not index.exists():
            raise FileNotFoundError("Provided FASTA index could not be found")

        # parse the index
        if isinstance(index, PathLike) and index.exists():
            with open(index) as index_file:
                self._parse_index(index_file)
        elif isinstance(index, TextIOBase):
            self._parse_index(index)
        else:
            self._index = None

        if index_compressed and isinstance(source, PathLike):
            # attempt to auto-locate the compressed index from the source path
            index_compressed = source.with_suffix(source.suffix + ".gzi")
        elif isinstance(index_compressed, PathLike) and not index_compressed.exists():
            raise FileNotFoundError("Provided compressed file index could not be found")
        elif isinstance(index_compressed, bool):
            index_compressed = None

        if index_compressed and (compression is not Compression.NONE):
            if compression is Compression.GZIP:
                raise ValueError(
                    "gzipped files can't have useful GZI indexes; use bgzip instead"
                )
            # parse the compressed index
            if isinstance(index_compressed, PathLike) and index_compressed.exists():
                with open(index_compressed, "rb") as index_compressed_file:
                    self._parse_index_compressed(index_compressed_file)
            elif isinstance(index_compressed, BufferedIOBase):
                self._parse_index_compressed(index_compressed)
        else:
            self._index_compressed = None

    def _parse_index(self, stream: TextIOBase):
        self._index = OrderedDict()

        for line in stream:
            fields = line.rstrip().split("\t")
            name = str(fields[0])
            self._index[name] = FaiEntry(
                name=name,
                length=int(fields[1]),
                offset=int(fields[2]),
                linebases=int(fields[3]),
                linewidth=int(fields[4]),
            )

    def _parse_index_compressed(self, stream: BufferedIOBase):
        blocks = [CompressedBlock(0, 0)]
        (number_entries,) = struct.unpack("<Q", stream.read(struct.calcsize("<Q")))
        for i in range(number_entries):
            compressed_offset, uncompressed_offset = struct.unpack(
                "<QQ", stream.read(struct.calcsize("<QQ"))
            )
            blocks.append(CompressedBlock(compressed_offset, uncompressed_offset))
        self._index_compressed = tuple(blocks)

    def _get_seqid_length(self, seqid: str) -> int:
        if self._index is not None:
            return self._index[seqid].length
        raise RuntimeError("Can't lazily get length of seqid without index")

    def __getitem__(self, key) -> FASTARecord:
        if self._index is None:
            raise RuntimeError("Random access to a FASTA file requires an index")

        if self._compression is Compression.GZIP:
            raise RuntimeError(
                "Random access into a plain gzipped FASTA file is not supported, use bgzip"
            )

        if isinstance(key, str):
            if key in self:
                return FASTARecord._from_FASTA_file(file=self, seqid=key)
            else:
                raise KeyError(f"SeqID '{key}' is not present in FASTA file")

        raise ValueError(
            f"FASTA files can only be indexed by strings, not '{type(key).__name__}'"
        )

    def __contains__(self, key) -> bool:
        if self._index is None:
            raise RuntimeError(
                "Checking if a sequence is in a FASTA file requires an index"
            )

        return key in self._index

    def close(self):
        self._stream.close()

    @property
    def closed(self) -> bool:
        return self._stream.closed

    def __enter__(self) -> "FASTAFile":
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _iter_indexed(self) -> Generator[FASTARecord, None, None]:
        seqids = list(self._index)
        for i, seqid in enumerate(seqids):
            # calculate start
            if i == 0:
                start = 0
            else:
                # compute total bytes past offset to start of next record
                previous = self._index[seqids[i - 1]]
                start = (
                    previous.offset
                    + (previous.linewidth * (previous.length // previous.linebases))
                    + (previous.length % previous.linebases)
                    + (previous.linewidth - previous.linebases)
                )

            # calculate bytes to read
            if (i + 1) == len(seqids):
                bytes_to_read = None
            else:
                index = self._index[seqid]
                bytes_to_read = (
                    (index.offset - start)
                    + (index.linewidth * (index.length // index.linebases))
                    + (index.length % index.linebases)
                )

            with self._stream_lock:
                self._stream.seek(start)
                record = self._stream.read(bytes_to_read)

            # parse record
            lines = record.split("\n")
            fields = lines[0].split(" ")
            seqid = fields[0][1:]
            description = " ".join(fields[1:]) if len(fields) > 1 else ""
            sequence = "".join(lines[1:])

            yield FASTARecord(id=seqid, sequence=sequence, description=description)

    def _iter_unindexed(self) -> Generator[FASTARecord, None, None]:
        current_offset = 0
        with self._stream_lock:
            self._stream.seek(current_offset)
            line = self._stream.readline()
            current_offset = self._stream.tell()

        if not line.startswith(">"):
            raise ValueError("First line in a FASTA file must start with '>'")

        while line:
            fields = line.rstrip()[1:].split(" ")
            seqid = fields[0]
            if len(fields) > 1:
                description = " ".join(fields[1:])
            else:
                description = ""

            sequence = []
            with self._stream_lock:
                self._stream.seek(current_offset)
                line = self._stream.readline()
                while line and (not line.startswith(">")):
                    sequence.append(line.rstrip())
                    line = self._stream.readline()
                current_offset = self._stream.tell()
            sequence = "".join(sequence)
            yield FASTARecord(id=seqid, sequence=sequence, description=description)

    def _iter_compressed(self) -> Generator[FASTARecord, None, None]:
        current_offset = 0
        with self._stream_lock:
            self._stream.seek(current_offset)
            buf = _read_chunk(self._stream, self._compression).decode("utf-8")
            current_offset = self._stream.tell()

        seqid = None
        while buf:
            if buf.startswith(">"):
                # yield current record, if any
                if seqid is not None:
                    yield FASTARecord(
                        id=seqid, sequence="".join(sequences), description=description
                    )

                description_line, _, buf = buf.partition("\n")
                fields = description_line[1:].split(" ")
                seqid = fields[0]
                description = " ".join(fields[1:]) if len(fields) > 1 else ""
                sequences = []
            else:
                sequence, sep, buf = buf.partition("\n")
                if sep:
                    sequences.append(sequence)
                else:
                    # no newline in block, append the rest
                    sequences.append(buf)
                    buf = ""

            if not buf:
                # maybe get next block
                with self._stream_lock:
                    self._stream.seek(current_offset)
                    buf += _read_chunk(self._stream, self._compression).decode("utf-8")
                    current_offset = self._stream.tell()

        if seqid is not None:
            yield FASTARecord(
                id=seqid, sequence="".join(sequences), description=description
            )

    def _iter_unindexed_gzip(self) -> Generator[FASTARecord, None, None]:
        with self._stream_lock:
            line = self._stream.readline()
            current_offset = self._stream.tell()

        seqid = None
        while line:
            if line.startswith(">"):
                if seqid is not None:
                    yield FASTARecord(
                        id=seqid, sequence="".join(sequences), description=description
                    )

                fields = line.rstrip()[1:].split(" ")
                seqid = fields[0]
                description = " ".join(fields[1:]) if len(fields) > 1 else ""
                sequences = []
            else:
                sequences.append(line.rstrip())

            with self._stream_lock:
                self._stream.seek(current_offset)
                line = self._stream.readline()
                current_offset = self._stream.tell()

        if seqid is not None:
            yield FASTARecord(
                id=seqid, sequence="".join(sequences), description=description
            )

    def __iter__(self) -> Generator[FASTARecord, None, None]:
        if self._compression is Compression.NONE:
            if self._index is None:
                return self._iter_unindexed()
            else:
                return self._iter_indexed()
        elif self._compression is Compression.GZIP:
            if self._index is None:
                return self._iter_unindexed_gzip()
            else:
                return self._iter_indexed()

        return self._iter_compressed()

    @property
    def closed(self):
        return self._stream.closed

    def _fetch(
        self,
        seqid: str,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
    ) -> str:
        if self.closed:
            raise ValueError("Cannot fetch a sequence from a closed FASTA file")

        if self._index is not None:
            index = self._index[seqid]
        else:
            raise RuntimeError("Cannot randomly access a FASTA file without an index")

        if (self._compression is not Compression.NONE) and (
            self._index_compressed is None
        ):
            raise RuntimeError(
                "Cannot randomly access a compressed FASTA file without a compressed index"
            )

        # set a default step of 1
        step = 1 if step is None else step

        if start is None:
            if step > 0:
                start = 0
            elif step < 0:
                start = index.length
        else:
            # convert negative to equivalent positive coordinate
            if start < 0:
                start = index.length + start
            # clamp start between (-1, len)
            start = max(-1, min(start, index.length))

        if stop is None:
            if step > 0:
                stop = index.length
            elif step < 0:
                stop = -1
        else:
            # convert negative to equivalent positive coordinate
            if stop < 0:
                stop = index.length + stop
            # clamp stop between (-1, len)
            stop = max(-1, min(stop, index.length))

        # handle zero-length requests
        if start == stop:
            return ""

        # handle "invalid" start/stop/step combinations
        if (step > 0) and (start > stop):
            return ""
        elif (step < 0) and (start < stop):
            return ""

        # step will be handled at the end
        if start > stop:
            start, stop = stop + 1, start + 1

        offset_start = self._base_to_byte_offset(
            start, index.linebases, index.linewidth
        )
        offset_stop = self._base_to_byte_offset(stop, index.linebases, index.linewidth)
        bytes_to_read = offset_stop - offset_start

        offset_start += index.offset
        if self._compression is Compression.NONE:
            with self._stream_lock:
                # seek to the start of the requested subsequence
                self._stream.seek(offset_start)
                sequence = "".join(self._stream.read(bytes_to_read).splitlines())
        else:
            # search for start block
            # BGZF has a defined maximum block size, can use that to limit search space
            min_start = (
                offset_start // _bgzf.MAX_BLOCK_SIZE_IN_BYTES
                if (self._compression is Compression.BGZF)
                else 0
            )
            block_start = bisect_left(
                a=self._index_compressed,
                x=offset_start,
                lo=min_start,
                key=attrgetter("uncompressed_offset"),
            ) - 1

            # scan for stop block
            block_stop = block_start + 1
            while (block_stop < len(self._index_compressed)) and (
                self._index_compressed[block_stop].uncompressed_offset
                < (offset_start + bytes_to_read)
            ):
                block_stop += 1

            uncompressed = ""
            compressed = []
            with self._stream_lock:
                self._stream.seek(self._index_compressed[block_start].compressed_offset)
                for i in range(block_start, block_stop):
                    if (i + 1) < len(self._index_compressed):
                        block_length = (
                            self._index_compressed[i + 1].compressed_offset
                            - self._index_compressed[i].compressed_offset
                        )
                        compressed.append(self._stream.read(block_length))
                    else:
                        # we've reached the last block, read the rest of the file
                        compressed.append(self._stream.read())

            uncompressed = _decompress(b"".join(compressed), self._compression).decode(
                "utf-8"
            )
            # subset the uncompressed data to the correct subsequence
            uncompressed_offset_start = (
                offset_start - self._index_compressed[block_start].uncompressed_offset
            )
            sequence = "".join(
                uncompressed[
                    uncompressed_offset_start : uncompressed_offset_start
                    + bytes_to_read
                ].split("\n")
            )

        return sequence[::step]

    @staticmethod
    def _base_to_byte_offset(
        base_offset: int, bases_per_line: int, bytes_per_line: int
    ) -> int:
        full_lines_before = base_offset // bases_per_line
        bases_before_on_line = base_offset % bases_per_line
        return (full_lines_before * bytes_per_line) + bases_before_on_line

    def __getstate__(self):
        if self._path is None:
            raise RuntimeError("Can't pickle a FASTAFile based on a stream")

        # https://docs.python.org/3/library/pickle.html#handling-stateful-objects
        state = self.__dict__.copy()

        del state["_stream"]
        del state["_stream_lock"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._stream = open(self._path)
        self._stream_lock = multiprocessing.Lock()
