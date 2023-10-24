import textwrap
import unittest
import zlib
from io import BytesIO, StringIO

import fastah.compression._bgzf as bgzf
import fastah.compression._zstd as _zstd
import fastah.fasta
from fastah import FASTAFile

if _zstd._HAS_ZSTD:
    import zstandard


class SimpleFASTA(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._source = textwrap.dedent(
            """\
            >seq1
            ACTG
            ACTG
            AC
            >seq2
            GTC
            G
            """
        )
        cls._index = textwrap.dedent(
            """\
            seq1	10	6	4	5
            seq2	4	25	3	4
            """
        )


class TestUncompressedFASTAIndexing(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    def test_index(self):
        index = StringIO("")
        fastah.fasta.index(StringIO(self._source), index)
        self.assertEqual(self._index, index.getvalue())


class TestUncompressedIndexedFASTAParsing(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls._file = FASTAFile(source=StringIO(cls._source), index=StringIO(cls._index))

    def test_index_int(self):
        self.assertEqual(self._file["seq1"][1], "C")
        self.assertEqual(self._file["seq2"][3], "G")
        self.assertEqual(self._file["seq1"][-2], "A")

    def test_index_int_invalid(self):
        with self.assertRaises(IndexError):
            self._file["seq1"][100]
        with self.assertRaises(IndexError):
            self._file["seq1"][10]
        with self.assertRaises(IndexError):
            self._file["seq1"][-11]
        with self.assertRaises(IndexError):
            self._file["seq2"][-5]

    def test_index_slice_single_line_small_range_first(self):
        self.assertEqual(self._file["seq1"][0:2], "AC")

    def test_index_slice_single_line_small_range_middle(self):
        self.assertEqual(self._file["seq1"][5:7], "CT")

    def test_index_slice_multi_line_end_at_start(self):
        self.assertEqual(self._file["seq1"][2:4], "TG")

    def test_index_slice_multi_line_end_in_middle(self):
        self.assertEqual(self._file["seq1"][2:6], "TGAC")

    def test_index_slice_long_stop(self):
        self.assertEqual(self._file["seq2"][:100], "GTCG")
        self.assertEqual(self._file["seq2"][:-10:-1], "GCTG")

    def test_index_slice_start_equal_stop(self):
        self.assertEqual(self._file["seq1"][4:4], "")

    def test_index_slice_start_greater_than_stop_but_positive_step(self):
        self.assertEqual(self._file["seq1"][6:4], "")

    def test_index_slice_start_less_than_stop_but_negative_step(self):
        self.assertEqual(self._file["seq1"][4:6:-1], "")

    def test_index_slice(self):
        with self.subTest(start=None, stop=None, step=None):
            self.assertEqual(self._file["seq1"][:], "ACTGACTGAC")
            self.assertEqual(self._file["seq2"][::], "GTCG")

        with self.subTest(start=None, stop=None, step="positive"):
            self.assertEqual(self._file["seq1"][::2], "ATATA")

        with self.subTest(start=None, stop=None, step="negative"):
            self.assertEqual(self._file["seq2"][::-1], "GCTG")

        with self.subTest(start=None, stop="positive", step=None):
            self.assertEqual(self._file["seq1"][:4], "ACTG")

        with self.subTest(start=None, stop="positive", step="positive"):
            self.assertEqual(self._file["seq1"][:4:2], "AT")

        with self.subTest(start=None, stop="positive", step="negative"):
            self.assertEqual(self._file["seq1"][:6:-2], "CG")

        with self.subTest(start=None, stop="negative", step=None):
            self.assertEqual(self._file["seq1"][:-2], "ACTGACTG")

        with self.subTest(start=None, stop="negative", step="positive"):
            self.assertEqual(self._file["seq1"][:-3:2], "ATAT")

        with self.subTest(start=None, stop="negative", step="negative"):
            self.assertEqual(self._file["seq1"][:-3:-1], "CA")

        with self.subTest(start="positive", stop=None, step=None):
            self.assertEqual(self._file["seq1"][8:], "AC")

        with self.subTest(start="positive", stop=None, step="positive"):
            self.assertEqual(self._file["seq1"][4::2], "ATA")

        with self.subTest(start="positive", stop=None, step="negative"):
            self.assertEqual(self._file["seq1"][7::-1], "GTCAGTCA")

        with self.subTest(start="positive", stop="positive", step=None):
            self.assertEqual(self._file["seq1"][6:8], "TG")

        with self.subTest(start="positive", stop="positive", step="positive"):
            self.assertEqual(self._file["seq1"][4:8:2], "AT")

        with self.subTest(start="positive", stop="positive", step="negative"):
            self.assertEqual(self._file["seq1"][8:4:-2], "AT")

        with self.subTest(start="positive", stop="negative", step=None):
            self.assertEqual(self._file["seq1"][7:-1], "GA")

        with self.subTest(start="positive", stop="negative", step="positive"):
            self.assertEqual(self._file["seq1"][6:-1:2], "TA")

        with self.subTest(start="positive", stop="negative", step="negative"):
            self.assertEqual(self._file["seq1"][8:-4:-1], "AG")

        with self.subTest(start="negative", stop=None, step=None):
            self.assertEqual(self._file["seq1"][-2:], "AC")

        with self.subTest(start="negative", stop=None, step="positive"):
            self.assertEqual(self._file["seq1"][-4::2], "TA")

        with self.subTest(start="negative", stop=None, step="negative"):
            self.assertEqual(self._file["seq1"][-6::-2], "ATA")

        with self.subTest(start="negative", stop="positive", step=None):
            self.assertEqual(self._file["seq1"][-4:8], "TG")

        with self.subTest(start="negative", stop="positive", step="positive"):
            self.assertEqual(self._file["seq1"][-4:8:2], "T")

        with self.subTest(start="negative", stop="positive", step="negative"):
            self.assertEqual(self._file["seq1"][-5:2:-2], "CG")

        with self.subTest(start="negative", stop="negative", step=None):
            self.assertEqual(self._file["seq1"][-5:-3], "CT")

        with self.subTest(start="negative", stop="negative", step="positive"):
            self.assertEqual(self._file["seq1"][-8:-1:2], "TATA")

        with self.subTest(start="negative", stop="negative", step="negative"):
            self.assertEqual(self._file["seq1"][-1:-8:-2], "CGCG")

    def test_index_slice_edge_cases(self):
        self.assertEqual(self._file["seq1"][-1:4], "")
        self.assertEqual(self._file["seq1"][-1:4:-1], "CAGTC")

    def test_index_slice_out_of_bounds(self):
        self.assertEqual(self._file["seq1"][15:20], "")
        self.assertEqual(self._file["seq1"][15:20:-1], "")
        self.assertEqual(self._file["seq1"][-20:-15], "")
        self.assertEqual(self._file["seq1"][-20:-15:-1], "")
        self.assertEqual(self._file["seq1"][-20:-15:-3], "")

    def test_record_length(self):
        self.assertEqual(len(self._file["seq1"]), 10)
        self.assertEqual(len(self._file["seq2"]), 4)

    def test_contains(self):
        self.assertTrue("seq1" in self._file)
        self.assertFalse("seqN" in self._file)

    def test_index_invalid_key(self):
        with self.assertRaises(KeyError):
            self._file["seqN"]
        with self.assertRaises(ValueError):
            self._file[42]


class TestBGZFCompressedIndexedFASTAParsing(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        compressed = BytesIO(b"")
        bgzf.compress(BytesIO(cls._source.encode("utf-8")), compressed)
        index_compressed = BytesIO(b"")
        bgzf.index(compressed, index_compressed)
        compressed.seek(0)
        index_compressed.seek(0)

        cls._file = FASTAFile(
            source=compressed,
            index=StringIO(cls._index),
            index_compressed=index_compressed,
        )

    def test_index_slice(self):
        self.assertEqual(self._file["seq1"][:4], "ACTG")


class TestUncompressedFASTAIterating(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls._file = FASTAFile(source=StringIO(cls._source))

    def test_iterating(self):
        records = list(self._file)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].sequence, "ACTGACTGAC")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].sequence, "GTCG")


class TestUncompressedIndexedFASTAIterating(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls._file = FASTAFile(source=StringIO(cls._source), index=StringIO(cls._index))

    def test_iterating(self):
        records = list(self._file)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].sequence, "ACTGACTGAC")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].sequence, "GTCG")


class TestBGZFCompressedIndexedFASTAIteration(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        compressed = BytesIO(b"")
        bgzf.compress(BytesIO(cls._source.encode("utf-8")), compressed)
        index_compressed = BytesIO(b"")
        bgzf.index(compressed, index_compressed)
        compressed.seek(0)
        index_compressed.seek(0)

        cls._file = FASTAFile(
            source=compressed,
            index=StringIO(cls._index),
            index_compressed=index_compressed,
        )

    def test_iterating(self):
        records = list(self._file)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].sequence, "ACTGACTGAC")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].sequence, "GTCG")


class TestGZIPCompressedFASTAIteration(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        compressed = BytesIO(zlib.compress(cls._source.encode("utf-8"), wbits=31))

        cls._file = FASTAFile(source=compressed)

    def test_iterating(self):
        records = list(self._file)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].sequence, "ACTGACTGAC")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].sequence, "GTCG")

    def test_index(self):
        with self.assertRaises(RuntimeError):
            self._file["seq1"][:4]


class TestGZIPCompressedFAIIndexedFASTAIteration(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        compressed = BytesIO(zlib.compress(cls._source.encode("utf-8"), wbits=31))

        cls._file = FASTAFile(source=compressed, index=StringIO(cls._index))

    def test_iterating(self):
        records = list(self._file)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].sequence, "ACTGACTGAC")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].sequence, "GTCG")

    def test_index(self):
        with self.assertRaises(RuntimeError):
            self._file["seq1"][:4]


@unittest.skipUnless(_zstd._HAS_ZSTD, "requires zstandard")
class TestZstdCompressedFAIIndexedFASTAIterating(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        compressor = zstandard.ZstdCompressor()
        compressed = BytesIO(compressor.compress(cls._source.encode("utf-8")))

        cls._file = FASTAFile(source=compressed, index=StringIO(cls._index))

    def test_iterating(self):
        records = list(self._file)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].sequence, "ACTGACTGAC")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].sequence, "GTCG")

    def test_index(self):
        with self.assertRaises(RuntimeError):
            self._file["seq1"][:4]


@unittest.skipUnless(_zstd._HAS_ZSTD, "requires zstandard")
class TestZStdCompressedIndexedFASTAParsing(SimpleFASTA):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        compressor = zstandard.ZstdCompressor()

        compressed = BytesIO(compressor.compress(cls._source.encode("utf-8")))
        index_compressed = BytesIO(b"")
        _zstd.index(compressed, index_compressed)
        compressed.seek(0)
        index_compressed.seek(0)

        cls._file = FASTAFile(
            source=compressed,
            index=StringIO(cls._index),
            index_compressed=index_compressed,
        )

    def test_index_slice(self):
        self.assertEqual(self._file["seq1"][:4], "ACTG")


if __name__ == "__main__":
    unittest.main()
