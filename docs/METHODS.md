# Методы и источники

Workflow обрабатывает публичные paired-end RRBS-данные: SRR17948848, SRR17948852, SRR17948854, SRR17948856. Источники FASTQ, контрольные суммы и группы образцов явно заданы в `config/samples.tsv`. Референсы и аннотации подготавливаются автоматически.

## Этапы обработки

| Практический этап | Реализация |
|---|---|
| FASTQ и sample sheet | Явный TSV с accession, группами, URL, MD5; автоматическая проверка |
| QC и trimming | FastQC до/после, Trim Galore paired `--rrbs --illumina`, MultiQC |
| Mapping | Bismark 0.25.1, Bowtie2 2.5.4, полный UCSC hg38 |
| Methylation extraction | CpG coverage, no-overlap, M-bias, Bismark HTML |
| Стандартные хромосомы | chr1–chr22, chrX, chrY; остальные исключаются только из downstream coverage |
| methylKit | mincov/filter/normalize/unite, гистограммы, Spearman, clustering, PCA, DMP |
| Региональный анализ | Независимый анализ фиксированных окон; не слияние DMP в de novo DMR |
| Геномная аннотация | Перекрытие UCSC генов и промоторов, CpG islands/shore 2 kb |
| Функциональная аннотация | GO ORA, измеренный фон, BH; метод отличается от ручного DAVID/g:Profiler |
| Интеграция с RNA | Дополнительный join с готовыми результатами RNA-seq, без обработки RNA FASTQ |

В PCA берутся до 100 000 вариабельных CpG с фиксированным seed. Дифференциальный анализ использует все подходящие позиции, а не подвыборку PCA. Удаление PC1 и клинические ковариаты по умолчанию не включены, поскольку их нельзя обосновать для этих четырёх записей.

Промотор: TSS −2000/+500 для плюсовой цепи, зеркально для минусовой; альтернативные транскрипты дают несколько возможных ассоциаций. Промотор имеет приоритет в объединённой колонке genes, отдельная колонка gene_body_genes сохраняет все пересечения с телом гена. CpG shore — 2 kb вокруг островков, исключая сами островки. Это координатная аннотация, не доказательство функционального эффекта.

GO анализирует уникальные SYMBOL, использует org.Hs.eg.db и GO.db из окружения. Термы включают ancestral GOALL. Фон — гены с измеренными позициями и конечным q-value, имеющие GO-аннотацию. Термы с 3–500 генами фона тестируются гипергеометрически, включая термы без foreground hits, затем BH по всем протестированным термам. Для окончательного исследования нужен метод, учитывающий специфическое для RRBS распределение числа CpG на ген.

UCSC FASTA и таблицы аннотаций доступны по изменяемым URL. При первом скачивании записываются SHA256, сами файлы нужно сохранять для точного повторения. Включён полный explicit lock проверенного Linux окружения; при изменении окружения сохраните новый lock. Повторное скачивание обновлённой аннотации может изменить результаты.

## Первичные источники

- [GEO GSE196387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196387) — исходное исследование.
- [ENA API: file reports](https://github.com/enasequence/read_docs/blob/master/retrieval/programmatic-access/file-reports.rst) — FASTQ, размеры и MD5.
- [Bismark 0.25.1](https://github.com/FelixKrueger/Bismark/tree/v0.25.1) — интерфейс и naming закреплённой версии.
- [Bismark: RRBS и дедупликация](https://felixkrueger.github.io/Bismark/usage/deduplication/) — почему coordinate dedup не применяется.
- [methylKit user guide](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html) — статистика и форматы.
- [methylKit calculateDiffMeth](https://github.com/al2na/methylKit/blob/master/man/calculateDiffMeth-methods.Rd) — overdispersion/test/adjust.
- [Snakemake deployment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) — окружения и повторяемость.
- [UCSC hg38 downloads](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/) — референс и аннотации.

