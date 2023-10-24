import argparse
import logging
from pathlib import Path

from . import compression
from .compression import Compression
from .fasta import FASTAFile
from .fasta import index as index_fasta


def index(args: argparse.Namespace):
    index = (
        args.path.with_suffix(args.path.suffix + ".fai")
        if args.index is None
        else args.index
    )
    if index.exists() and not args.force:
        raise FileExistsError(
            "FASTA index already exists; use --force if you really want to overwrite it"
        )

    index_fasta(args.path, index)

    index_compressed = (
        args.path.with_suffix(args.path.suffix + ".gzi")
        if args.index_compressed is None
        else args.index_compressed
    )

    # if compressed, create GZI index
    with open(args.path, "rb") as FASTA:
        try:
            compression_format = compression._detect_format(FASTA)
            if compression_format is Compression.GZIP:
                logging.warning(
                    "GZI indices can't be created for plain gzip files, skipping"
                )
            else:
                with open(index_compressed, "wb") as index_compressed_file:
                    compression._index(FASTA, index_compressed_file, compression_format)
        except RuntimeError:
            # not compressed in any format we know how to handle, ignore
            pass


def _parse_region(region: str) -> tuple[str, int, int]:
    seqid, startstop = region.split(":")
    start, stop = startstop.split("-")

    return seqid, int(start) - 1, int(stop)


def fetch(args: argparse.Namespace):
    seqid, start, stop = _parse_region(args.region)

    index = (
        args.path.with_suffix(args.path.suffix + ".fai")
        if args.index is None
        else args.index
    )

    index_compressed = (
        args.path.with_suffix(args.path.suffix + ".gzi")
        if args.index_compressed is None
        else args.index_compressed
    )

    with FASTAFile(
        source=args.path, index=index, index_compressed=index_compressed
    ) as reference:
        print(reference[seqid][start:stop])


def main():
    parser = argparse.ArgumentParser(
        description="A toolkit for operating on FASTA sequence files"
    )
    subparsers = parser.add_subparsers(
        title="subcommands", metavar="subcommand", required=True
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="output progress and other informative messages",
        action="count",
        dest="verbosity",
        default=0,
    )

    parser_index = subparsers.add_parser("index", help="index FASTA files")
    parser_index.add_argument(
        "path", help="path to the FASTA file to be indexed", type=Path
    )
    parser_index.add_argument(
        "-i", "--index", help="path to write the FASTA index", type=Path
    )
    parser_index.add_argument(
        "-z",
        "--index_compressed",
        help="path to write the compressed file index",
        type=Path,
    )
    parser_index.add_argument(
        "-f",
        "--force",
        help="overwrite indicies if they already exist",
        action="store_true",
    )
    parser_index.set_defaults(target=index)

    parser_fetch = subparsers.add_parser(
        "fetch", help="fetch sequences from FASTA files"
    )
    parser_fetch.add_argument(
        "path", help="path to the FASTA file to be queried", type=Path
    )
    parser_fetch.add_argument(
        "region",
        help="region to extract from FASTA file in samtools format (seqid:start-end, 1-based closed intervals)",
    )
    parser_fetch.add_argument(
        "-i", "--index", help="path to the FASTA index (.fai)", type=Path
    )
    parser_fetch.add_argument(
        "-z",
        "--index_compressed",
        help="path to the compressed file index (.gzi)",
        type=Path,
    )
    parser_fetch.set_defaults(target=fetch)

    args = parser.parse_args()

    logging_level = max(logging.DEBUG, logging.WARNING - (args.verbosity * 10))
    logging.basicConfig(
        format="{asctime} [{module}:{levelname}] {message}"
        if logging_level is logging.DEBUG
        else "{asctime} {levelname}: {message}",
        style="{",
        level=logging_level,
        force=True,
    )

    args.target(args)
