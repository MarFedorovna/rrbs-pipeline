suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(jsonlite))
cfg <- snakemake@params[["cfg"]]
out <- snakemake@params[["out"]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
logcon <- file(snakemake@log[[1]], "wt")
sink(logcon); sink(logcon, type="message")
bed <- function(path) {
    d <- read.delim(path, header=FALSE, stringsAsFactors=FALSE)
    makeGRangesFromDataFrame(data.frame(chr=d[[1]], start=d[[2]]+1L, end=d[[3]],
                                       feature=d[[4]]), keep.extra.columns=TRUE)
}
genes <- bed(snakemake@input[["genes"]])
promoters <- bed(snakemake@input[["promoters"]])
islands <- bed(snakemake@input[["islands"]])
shores <- resize(islands, width=width(islands)+4000L, fix="center")
start(shores) <- pmax(1L, start(shores))
overlap_names <- function(query, subject) {
    hits <- findOverlaps(query, subject, ignore.strand=TRUE)
    result <- rep("", length(query))
    if (length(hits)) {
        labels <- tapply(subject$feature[subjectHits(hits)], queryHits(hits),
                          function(x) paste(sort(unique(x)), collapse=";"))
        result[as.integer(names(labels))] <- unname(labels)
    }
    result
}
split_genes <- function(values) unique(unlist(strsplit(values[nzchar(values)], ";", fixed=TRUE)))
go_test <- function(selected_symbols, background_symbols) {
    empty <- data.frame(GO=character(), ontology=character(), term=character(), selected_hits=integer(),
                        background_hits=integer(), selected_total=integer(), background_total=integer(),
                        pvalue=numeric(), padj=numeric())
    if (!isTRUE(cfg$go_enrichment)) return(empty)
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    suppressPackageStartupMessages(library(GO.db))
    suppressPackageStartupMessages(library(AnnotationDbi))
    background_symbols <- intersect(background_symbols, keys(org.Hs.eg.db, keytype="SYMBOL"))
    selected_symbols <- intersect(selected_symbols, background_symbols)
    if (length(selected_symbols) < 2 || length(background_symbols) < 3) return(empty)
    map <- suppressMessages(select(org.Hs.eg.db, keys=background_symbols,
                                   keytype="SYMBOL", columns=c("GOALL", "ONTOLOGYALL")))
    map <- unique(map[!is.na(map$GOALL), c("SYMBOL", "GOALL", "ONTOLOGYALL")])
    # Restrict the universe to measured genes with an available GO annotation.
    background_symbols <- intersect(background_symbols, map$SYMBOL)
    selected_symbols <- intersect(selected_symbols, background_symbols)
    if (length(selected_symbols) < 2) return(empty)
    sets <- split(map$SYMBOL, map$GOALL)
    sets <- lapply(sets, unique)
    sets <- sets[lengths(sets) >= cfg$go_min_genes & lengths(sets) <= cfg$go_max_genes]
    if (!length(sets)) return(empty)
    n <- length(selected_symbols); N <- length(background_symbols)
    k <- vapply(sets, function(x) sum(x %in% selected_symbols), integer(1))
    m <- lengths(sets)
    result <- data.frame(GO=names(sets), selected_hits=k, background_hits=m,
                         selected_total=n, background_total=N,
                         pvalue=phyper(k-1L, m, N-m, n, lower.tail=FALSE))
    result$padj <- p.adjust(result$pvalue, "BH")
    descriptions <- suppressMessages(select(GO.db, keys=result$GO, columns=c("TERM", "ONTOLOGY"), keytype="GOID"))
    result$term <- descriptions$TERM[match(result$GO, descriptions$GOID)]
    result$ontology <- descriptions$ONTOLOGY[match(result$GO, descriptions$GOID)]
    result[order(result$padj, result$pvalue), ]
}
summaries <- list()
for (prefix in c("dmp", "tiles")) {
    path <- snakemake@input[[if (prefix == "dmp") "dmp" else "tiles"]]
    table <- read.delim(gzfile(path), stringsAsFactors=FALSE)
    query <- makeGRangesFromDataFrame(table, keep.extra.columns=FALSE)
    table$promoter_genes <- overlap_names(query, promoters)
    table$gene_body_genes <- overlap_names(query, genes)
    table$cpg_context <- ifelse(countOverlaps(query, islands, ignore.strand=TRUE) > 0, "island",
        ifelse(countOverlaps(query, shores, ignore.strand=TRUE) > 0, "shore_2kb", "other"))
    # Promoter associations have priority; gene body overlaps remain separately available.
    table$genes <- ifelse(nzchar(table$promoter_genes), table$promoter_genes, table$gene_body_genes)
    table$significant <- is.finite(table$qvalue) & table$qvalue < cfg$qvalue &
        is.finite(table$meth.diff) & abs(table$meth.diff) > cfg$difference
    con <- gzfile(file.path(out, paste0(prefix, "_annotated.tsv.gz")), "wt")
    write.table(table, con, sep="\t", quote=FALSE, row.names=FALSE); close(con)
    foreground <- split_genes(table$genes[table$significant])
    background <- split_genes(table$genes[is.finite(table$qvalue)])
    writeLines(sort(foreground), file.path(out, paste0(prefix, "_genes.txt")))
    writeLines(sort(background), file.path(out, paste0(prefix, "_background_genes.txt")))
    enrichment <- go_test(foreground, background)
    write.table(enrichment, file.path(out, paste0(prefix, "_go.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
    summaries[[prefix]] <- list(tested=nrow(table), significant=sum(table$significant),
        foreground_genes=length(foreground), background_genes=length(background),
        go_terms_tested=nrow(enrichment), go_terms_fdr05=sum(enrichment$padj < 0.05, na.rm=TRUE))
}
write_json(list(results=summaries,
    annotation="UCSC ncbiRefSeqCurated/hg38; overlapping promoters and gene bodies; CpG island and 2kb shores",
    go="Exploratory hypergeometric ORA with BH across all eligible terms, using tested annotated genes as universe. RRBS CpG density and gene length bias are not corrected.",
    note="An overlapping gene is an association, not evidence of regulation or causality."),
    file.path(out, "summary.json"), pretty=TRUE, auto_unbox=TRUE)
sink(type="message"); sink(); close(logcon)
