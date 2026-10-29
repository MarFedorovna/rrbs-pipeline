suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(jsonlite))
cfg <- snakemake@params[["cfg"]]
out <- snakemake@params[["out"]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
logcon <- file(snakemake@log[[1]], "wt")
sink(logcon); sink(logcon, type="message")
set.seed(42)
meta <- read.delim(snakemake@input[["samples"]], stringsAsFactors=FALSE)
meta <- meta[match(unlist(snakemake@params[["ids"]]), meta$sample), ]
stopifnot(!anyNA(meta$sample), all(table(meta$treatment) >= 2))
write.table(meta[, c("sample", "run", "group", "treatment")],
            file.path(out, "sample_order.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

raw <- methRead(as.list(unlist(snakemake@input[["coverage"]])),
                sample.id=as.list(meta$sample), assembly=cfg$assembly,
                header=FALSE, pipeline="bismarkCoverage", context="CpG",
                treatment=as.integer(meta$treatment), mincov=cfg$min_coverage)
filtered <- filterByCoverage(raw, lo.count=cfg$min_coverage,
    hi.perc=if (cfg$max_coverage_percentile < 100) cfg$max_coverage_percentile else NULL)
counts <- data.frame(sample=meta$sample,
    loaded=vapply(raw, function(x) nrow(getData(x)), integer(1)),
    filtered=vapply(filtered, function(x) nrow(getData(x)), integer(1)))
write.table(counts, file.path(out, "site_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (any(counts$filtered == 0)) stop("Coverage filtering removed every CpG in a sample")

pdf(file.path(out, "qc.pdf"), width=8, height=6)
for (i in seq_along(raw)) {
    getMethylationStats(raw[[i]], plot=TRUE, both.strands=FALSE)
    title(main=paste(meta$sample[i], "methylation"))
    getCoverageStats(raw[[i]], plot=TRUE, both.strands=FALSE)
    title(main=paste(meta$sample[i], "coverage"))
}
normalized <- if (isTRUE(cfg$normalize_coverage)) normalizeCoverage(filtered) else filtered
meth <- unite(normalized, destrand=FALSE, min.per.group=as.integer(cfg$min_per_group))
if (nrow(getData(meth)) < 3) stop("Fewer than three shared CpGs remain; check data and groups")
saveRDS(meth, file.path(out, "methylBase.rds"))
matrix <- percMethylation(meth)
# Use complete CpGs for exploratory QC when a larger cohort permits missing values.
matrix <- matrix[complete.cases(matrix), , drop=FALSE]
if (nrow(matrix) < 3) stop("Not enough complete CpGs for QC plots")
colnames(matrix) <- meta$sample
correlation <- cor(matrix, method="spearman")
write.table(correlation, file.path(out, "correlation.tsv"), sep="\t", quote=FALSE)
heatmap(correlation, symm=TRUE, main="Spearman sample correlation")
plot(hclust(as.dist(1-correlation), method="average"), main="Sample clustering (1 - Spearman)")
varying <- which(apply(matrix, 1, var) > 0)
if (length(varying) >= 2) {
    selected <- if (length(varying) > 100000) sample(varying, 100000) else varying
    pca <- prcomp(t(matrix[selected, , drop=FALSE]), center=TRUE, scale.=FALSE)
    scores <- data.frame(sample=meta$sample, group=meta$group, pca$x)
    write.table(scores, file.path(out, "pca_scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    plot(pca$x[, 1], pca$x[, 2], col=as.integer(meta$treatment)+1, pch=19,
         xlab=paste0("PC1 (", variance[1], "%)"), ylab=paste0("PC2 (", variance[2], "%)"),
         main="PCA (up to 100,000 variable CpGs; seed 42)")
    text(pca$x[, 1], pca$x[, 2], labels=meta$sample, pos=3, cex=0.65)
} else {
    write.table(data.frame(sample=meta$sample, group=meta$group, PC1=NA, PC2=NA),
                file.path(out, "pca_scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    plot.new(); text(0.5, 0.5, "PCA unavailable: fewer than two variable CpGs")
}
dev.off()
rm(matrix, raw, filtered); gc()

write_gz <- function(table, filename) {
    con <- gzfile(file.path(out, filename), "wt")
    write.table(table, con, sep="\t", row.names=FALSE, quote=FALSE)
    close(con)
}
empty <- data.frame(chr=character(), start=integer(), end=integer(), strand=character(),
                    pvalue=numeric(), qvalue=numeric(), meth.diff=numeric())
analyze <- function(object, prefix) {
    if (is.null(object) || nrow(getData(object)) < 1) {
        result <- empty
    } else {
        difference <- calculateDiffMeth(object, adjust=cfg$adjust,
            overdispersion=cfg$overdispersion, test=cfg$test,
            mc.cores=as.integer(snakemake@threads), effect="wmean")
        result <- getData(difference)
    }
    keep <- is.finite(result$qvalue) & result$qvalue < cfg$qvalue &
            is.finite(result$meth.diff) & abs(result$meth.diff) > cfg$difference
    significant <- result[keep, , drop=FALSE]
    write_gz(result, paste0(prefix, "_all.tsv.gz"))
    write_gz(significant, paste0(prefix, "_significant.tsv.gz"))
    write_gz(significant[significant$meth.diff > 0, , drop=FALSE], paste0(prefix, "_hyper.tsv.gz"))
    write_gz(significant[significant$meth.diff < 0, , drop=FALSE], paste0(prefix, "_hypo.tsv.gz"))
    # TSV coordinates are one-based inclusive. Exported BED is zero-based half-open.
    bed <- data.frame(chr=significant$chr, start=significant$start-1L, end=significant$end,
                      id=if (nrow(significant)) paste0(prefix, "_", seq_len(nrow(significant))) else character())
    write.table(bed, file.path(out, paste0(prefix, "_significant.bed")),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    list(tested=nrow(result), significant=nrow(significant),
         hyper=sum(significant$meth.diff > 0), hypo=sum(significant$meth.diff < 0),
         nonfinite_q=sum(!is.finite(result$qvalue)))
}
base_summary <- analyze(meth, "dmp")
tile_summary <- if (isTRUE(cfg$tiles)) {
    tiles <- tileMethylCounts(normalized, win.size=as.integer(cfg$tile_size),
        step.size=as.integer(cfg$tile_step), cov.bases=as.integer(cfg$tile_min_cpg))
    united_tiles <- unite(tiles, destrand=FALSE, min.per.group=as.integer(cfg$min_per_group))
    analyze(united_tiles, "tiles")
} else analyze(NULL, "tiles")
summary <- list(assembly=cfg$assembly, samples=meta$sample,
    comparison=paste(unique(meta$group[meta$treatment == 1]), "minus", unique(meta$group[meta$treatment == 0])),
    effect_units="percentage points", qvalue_cutoff=cfg$qvalue, absolute_difference_cutoff=cfg$difference,
    test=cfg$test, overdispersion=cfg$overdispersion, adjust=cfg$adjust,
    dmp=base_summary, tiles=tile_summary,
    note="Two replicates per group: exploratory inference. Significant tiles are fixed windows, not de novo DMRs.")
write_json(summary, file.path(out, "summary.json"), pretty=TRUE, auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
sink(type="message"); sink(); close(logcon)
