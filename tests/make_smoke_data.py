#!/usr/bin/env python3
"""Small, artificial RRBS-like paired reads. Never biological evidence."""
import csv
import gzip
import hashlib
import random
from pathlib import Path

import yaml

root = Path("tests/generated")
root.mkdir(parents=True, exist_ok=True)
rng = random.Random(42)
fragments = ["CCGG" + "".join(rng.choices("ACGT", k=242)) + "CCGG" for _ in range(50)]
genome = "".join(fragments)
with gzip.open(root / "genome.fa.gz", "wt") as handle:
    handle.write(">chr1\n" + genome + "\n")
rows = []
complement = str.maketrans("ACGT", "TGCA")
for sid, treatment, probability in [("test1", 1, .85), ("test2", 1, .8), ("ctrl1", 0, .2), ("ctrl2", 0, .15)]:
    paths = [root / f"{sid}_R{mate}.fastq.gz" for mate in [1, 2]]
    with gzip.open(paths[0], "wt") as r1, gzip.open(paths[1], "wt") as r2:
        for index, fragment in enumerate(fragments):
            for copy in range(40):
                # Both mates come from one bisulfite-converted top-strand molecule.
                converted = "".join("C" if base == "C" and i + 1 < len(fragment) and fragment[i+1] == "G" and rng.random() < probability
                                    else "T" if base == "C" else base
                                    for i, base in enumerate(fragment))
                # MspI cuts C|CGG. Include the RRBS-compatible CGG start.
                left = converted[1:101]
                right = converted[-101:-1].translate(complement)[::-1]
                readid = f"{sid}:{index}:{copy}"
                r1.write(f"@{readid}/1\n{left}\n+\n{'I'*len(left)}\n")
                r2.write(f"@{readid}/2\n{right}\n+\n{'I'*len(right)}\n")
    rows.append(dict(sample=sid, run="synthetic", group="test" if treatment else "control", treatment=treatment,
                     r1=str(paths[0]), r2=str(paths[1]),
                     md5_r1=hashlib.md5(paths[0].read_bytes()).hexdigest(),
                     md5_r2=hashlib.md5(paths[1].read_bytes()).hexdigest()))
with open(root / "samples.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
with gzip.open(root / "genes.txt.gz", "wt") as handle:
    for i, symbol in enumerate(["TP53", "EGFR", "BRCA1", "BRCA2", "MYC"]):
        start, end = i * 2500, (i + 1) * 2500
        record = [0, "synthetic_tx_"+str(i), "chr1", "+", start, end, start, end, 1,
                  str(start)+",", str(end)+",", 0, symbol, "cmpl", "cmpl", "0,"]
        handle.write("\t".join(map(str, record)) + "\n")
with gzip.open(root / "islands.txt.gz", "wt") as handle:
    handle.write("0\tchr1\t100\t1000\tCpG:synthetic\t900\t80\t200\t20\t40\t1\n")
cfg = yaml.safe_load(Path("config/config.yaml").read_text())
cfg.update(samples=str(root / "samples.tsv"), outdir="results/smoke", reference_dir="reference/smoke",
           reference_url=str(root / "genome.fa.gz"), genes_url=str(root / "genes.txt.gz"),
           islands_url=str(root / "islands.txt.gz"), adjust="BH", min_coverage=5,
           max_coverage_percentile=100, tile_min_cpg=2, canonical_chromosomes=["chr1"],
           go_enrichment=True)
cfg["threads"] = dict(index=2, trim=4, align=4, extract=3, differential=2)
cfg["memory_mb"] = dict(index=1000, trim=2000, align=2000, extract=3000, differential=4000)
(root / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
print("Synthetic test fixtures generated. No human reads were downloaded.")
