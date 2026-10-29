#!/usr/bin/env bash
set -euo pipefail
printf 'OS / CPU / RAM / disk\n'
uname -a
lscpu | head -20
free -h
df -h .
printf '\nRequired tools\n'
for tool in snakemake python Rscript fastqc trim_galore bismark bismark_genome_preparation bismark_methylation_extractor bismark2report bowtie2 samtools multiqc bedtools; do
  command -v "$tool" || exit 1
done
Rscript -e 'stopifnot(requireNamespace("methylKit"), requireNamespace("org.Hs.eg.db"), requireNamespace("GO.db")); sessionInfo()'
printf '\nGPU is not used by this pipeline.\n'
