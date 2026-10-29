#!/usr/bin/env python3
"""Join promoter-associated DMPs to an existing RNA-seq differential expression TSV."""
import argparse
import csv
import gzip
from collections import defaultdict
from statistics import mean

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("annotated_dmp", help="annotation/dmp_annotated.tsv.gz")
p.add_argument("rna_tsv", help="One row per gene; columns gene, log2FoldChange, padj")
p.add_argument("output_tsv")
a = p.parse_args()
genes = defaultdict(list)
with gzip.open(a.annotated_dmp, "rt") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        for gene in filter(None, row["promoter_genes"].split(";")):
            genes[gene].append(row)
with open(a.rna_tsv) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if not {"gene", "log2FoldChange", "padj"}.issubset(reader.fieldnames or []):
        raise ValueError("RNA TSV needs gene, log2FoldChange, padj")
    rna = {}
    for row in reader:
        if row["gene"] in rna:
            raise ValueError("RNA table must have unique gene symbols; aggregate transcripts first")
        rna[row["gene"]] = row
with open(a.output_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene", "promoter_tested_cpg", "promoter_significant_cpg",
                     "mean_tested_cpg_methylation_difference_pp", "rna_log2FoldChange", "rna_padj"])
    for gene in sorted(genes.keys() & rna.keys()):
        rows = genes[gene]
        writer.writerow([gene, len(rows), sum(r["significant"] == "TRUE" for r in rows),
                         mean(float(r["meth.diff"]) for r in rows),
                         rna[gene]["log2FoldChange"], rna[gene]["padj"]])
