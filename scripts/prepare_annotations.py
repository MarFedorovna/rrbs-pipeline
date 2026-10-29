#!/usr/bin/env python3
"""UCSC genePred and cpgIslandExt to BED; coordinates stay zero-based half-open."""
import argparse
import gzip
from pathlib import Path


def prepare(genes_path, islands_path, outdir, allowed, upstream, downstream):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    genes, promoters, islands = set(), set(), set()
    with gzip.open(genes_path, "rt") as handle:
        for line in handle:
            a = line.rstrip().split("\t")
            chrom, strand, start, end, symbol = a[2], a[3], int(a[4]), int(a[5]), a[12]
            if chrom not in allowed:
                continue
            genes.add((chrom, start, end, symbol, ".", strand))
            tss = start if strand == "+" else end
            left, right = ((tss - upstream, tss + downstream) if strand == "+"
                           else (tss - downstream, tss + upstream))
            promoters.add((chrom, max(0, left), right, symbol, ".", strand))
    with gzip.open(islands_path, "rt") as handle:
        for line in handle:
            a = line.rstrip().split("\t")
            if a[1] in allowed:
                islands.add((a[1], int(a[2]), int(a[3]), a[4]))
    for name, rows in [("genes", genes), ("promoters", promoters), ("islands", islands)]:
        if not rows:
            raise ValueError(f"No {name} features match the selected chromosome names")
        with open(out / (name + ".bed"), "w") as handle:
            for row in sorted(rows):
                handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genes"); p.add_argument("islands"); p.add_argument("outdir")
    p.add_argument("--chromosomes", nargs="+", required=True)
    p.add_argument("--upstream", type=int, default=2000)
    p.add_argument("--downstream", type=int, default=500)
    a = p.parse_args()
    prepare(a.genes, a.islands, a.outdir, set(a.chromosomes), a.upstream, a.downstream)
