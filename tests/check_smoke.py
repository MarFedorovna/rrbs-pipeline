import gzip
import json
from pathlib import Path

root = Path("results/smoke")
summary = json.loads((root / "analysis/summary.json").read_text())
assert summary["dmp"]["tested"] > 50, summary
assert summary["dmp"]["hyper"] > 0, summary
assert summary["dmp"]["hyper"] > summary["dmp"]["hypo"], summary
assert summary["tiles"]["tested"] > 0, summary
assert (root / "multiqc/multiqc_report.html").stat().st_size > 1000
assert (root / "analysis/qc.pdf").read_bytes().startswith(b"%PDF")
assert len(list((root / "bismark").glob("*/report.html"))) == 4
annotation = json.loads((root / "annotation/summary.json").read_text())
assert annotation["results"]["dmp"]["background_genes"] > 0
assert annotation["results"]["dmp"]["go_terms_tested"] > 0
with gzip.open(root / "annotation/dmp_annotated.tsv.gz", "rt") as handle:
    assert "promoter_genes" in handle.readline()
provenance = json.loads((root / "provenance/run.json").read_text())
assert "0.25.1" in provenance["software"]["bismark"]
print("Real-tool synthetic end-to-end smoke checks passed.")
