#!/usr/bin/env python3
import argparse
import gzip


def filter_coverage(source, target, allowed):
    count = 0
    with gzip.open(source, "rt") as inp, gzip.open(target, "wt") as out:
        for line in inp:
            fields = line.rstrip().split("\t")
            if fields[0] not in allowed:
                continue
            if len(fields) != 6 or int(fields[1]) < 1 or int(fields[1]) != int(fields[2]):
                raise ValueError("Expected six-column, one-based Bismark CpG coverage")
            out.write(line)
            count += 1
    if not count:
        raise ValueError("No CpG coverage rows on allowed chromosomes; check reference naming and alignment")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("source"); p.add_argument("target"); p.add_argument("chromosomes", nargs="+")
    a = p.parse_args()
    filter_coverage(a.source, a.target, set(a.chromosomes))
