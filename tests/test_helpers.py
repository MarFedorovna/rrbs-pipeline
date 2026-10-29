import csv
import gzip
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def module(path):
    spec = importlib.util.spec_from_file_location(Path(path).stem, ROOT / path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


lib = module("workflow/lib.py")
fetch = module("scripts/fetch.py")
coverage = module("scripts/filter_coverage.py")
annotations = module("scripts/prepare_annotations.py")


class DataSafetyTests(unittest.TestCase):
    def test_sample_groups_and_order(self):
        samples = lib.load_samples(ROOT / "config/samples.tsv")
        self.assertEqual([s["run"] for s in samples.values()],
                         ["SRR17948848", "SRR17948852", "SRR17948854", "SRR17948856"])
        self.assertEqual([s["treatment"] for s in samples.values()], ["1", "0", "0", "1"])
        self.assertEqual({s["group"] for s in samples.values()}, {"Active", "T2D"})

    def test_duplicate_id_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "samples.tsv"
            p.write_text((ROOT / "config/samples.tsv").read_text().replace("Active_RES09", "T2D_RES20"))
            with self.assertRaisesRegex(ValueError, "duplicate"):
                lib.load_samples(p)

    def test_bad_checksum_cannot_publish_output(self):
        with tempfile.TemporaryDirectory() as d:
            source, target = Path(d) / "in", Path(d) / "out"
            source.write_bytes(b"abc")
            original_sleep = fetch.time.sleep
            fetch.time.sleep = lambda _: None
            try:
                with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
                    fetch.fetch(str(source), target, "0"*32)
                self.assertFalse(target.exists())
                self.assertFalse(target.with_name("out.part").exists())
                fetch.fetch(str(source), target, hashlib.md5(b"abc").hexdigest())
                self.assertEqual(target.read_bytes(), b"abc")
                source.write_bytes(b"changed")
                with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
                    fetch.fetch(str(source), target, hashlib.md5(b"abc").hexdigest())
            finally:
                fetch.time.sleep = original_sleep

    def test_chromosome_filter_and_one_based_positions(self):
        with tempfile.TemporaryDirectory() as d:
            source, target = Path(d) / "in.gz", Path(d) / "out.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("chr1\t1\t1\t50\t5\t5\nchr1_alt\t2\t2\t0\t0\t10\n")
            coverage.filter_coverage(source, target, {"chr1"})
            with gzip.open(target, "rt") as handle:
                self.assertEqual(handle.read(), "chr1\t1\t1\t50\t5\t5\n")

    def test_promoters_respect_both_strands(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            with gzip.open(root / "genes.gz", "wt") as handle:
                for strand in ["+", "-"]:
                    handle.write("\t".join(map(str, [0,"tx","chr1",strand,100,500,100,500,1,"100,","500,",0,strand]))+"\n")
            with gzip.open(root / "islands.gz", "wt") as handle:
                handle.write("0\tchr1\t100\t200\tisland\n")
            annotations.prepare(root / "genes.gz", root / "islands.gz", root / "out", {"chr1"}, 200, 50)
            rows = list(csv.reader((root / "out/promoters.bed").read_text().splitlines(), delimiter="\t"))
            self.assertEqual({r[3]: (int(r[1]), int(r[2])) for r in rows}, {"+": (0,150), "-": (450,700)})


if __name__ == "__main__":
    unittest.main()
