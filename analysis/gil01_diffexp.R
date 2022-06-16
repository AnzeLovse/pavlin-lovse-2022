library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(RUVSeq)
library(RColorBrewer)
library(EnhancedVolcano)

source("analysis/helper_functions.R")

RAW_DATA <- file.path("data", "raw")
BASE <- file.path("results", "gil01_experiment")
dir.create(BASE, showWarnings = FALSE)

counts <- read_counts(path = file.path(RAW_DATA, "gil01_experiment", "merged_counts.tsv"))
coldata <- read_metadata(
  samples = file.path(RAW_DATA, "gil01_experiment", "samples.tsv"),
  units = file.path(RAW_DATA, "gil01_experiment", "units.tsv")
)

# Make the column order match.
coldata <- coldata[colnames(counts), ]

# Filter genes by requiring more than 5 reads in at least two samples.
filtered <- counts[rowSums(counts > 5) >= 2, ]

# Plot RLE.
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = coldata)
set <- betweenLaneNormalization(set, which = "upper")
colors <- brewer.pal(6, "Set2")

pdf(file.path(BASE, "RLE.pdf"), width = 12, height = 8)
plotRLE(
  x = set,
  outline = FALSE,
  ylim = c(-4, 4),
  col = colors[coldata$condition]
)
dev.off()

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = coldata, design = ~ 0 + condition)
dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE) # variance stabilising transformation
pca_data <- plotPCA(vsd, intgroup = c("condition"), ntop = nrow(filtered), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
ggsave(file.path(BASE, "pca.pdf"))

comparisons <- list(
  G_M_vs_B_M = c("G_M", "B_M"),
  G_vs_B = c("G", "B"),
  B_M_vs_B = c("B_M", "B")
)

labels <- list(
  G_M_vs_B_M = c(""),
  G_vs_B = labels <- c(""),
  B_M_vs_B = c(
    "HIS92_RS04150", "HIS92_RS04600", "HIS92_RS06675", "HIS92_RS11675",
    "HIS92_RS11680", "HIS92_RS03705", "HIS92_RS04025", "HIS92_RS15270"
  )
)

result_list <- list()

for (name in names(comparisons)) {
  result <- results(dds, contrast = c("condition", comparisons[[name]]), alpha = 0.05)
  result <- result[order(result$padj), ]

  save_deg(
    result = result,
    out_file = file.path(BASE, paste0(name, ".tsv")),
    out_annotated = file.path(BASE, paste0(name, "_annotated.tsv")),
    fdr_threshold = 0.05
  )

  result_list[[name]] <- result
}

# Remove gil genes for plotting purposes.
gene_annotation <- read.csv(
  file = file.path("data/external/annotation_with_gil01.bed"),
  sep = "\t",
  col.names = c("chr", "start", "end", "gene", "type", "strand"),
  header = FALSE
)
gil_ids <- gene_annotation[gene_annotation$chr == "AJ536073.2", "gene"]

filtered_results <- lapply(result_list, function(x) {
  x[!rownames(x) %in% gil_ids, ]
})


xmax <- max(unlist(lapply(filtered_results, function(x) {
  max(abs(x$log2FoldChange))
})))
xlim <- c(-xmax, xmax)

gene_metadata <- get_gene_annot()

for (name in names(result_list)) {
  volcano_title <- stringr::str_replace(string = name, pattern = "_", replacement = " ")
  plot_volcano(
    result = filtered_results[[name]],
    gene_metadata = gene_metadata,
    outfile = file.path(BASE, paste0(name, "_volcano.pdf")),
    xlim = xlim,
    title = volcano_title,
    fdr_threshold = 0.05,
    labels = labels[[name]]
  )
}

# Approximate the amount of reads on the regulatory region of GIL01.
norm_counts <- counts(dds, normalized = TRUE)
gil_rows <- rownames(norm_counts) %in% gil_ids[1:9]
gm_gil_reads <- sum(norm_counts[gil_rows, c("G1M.t0.1", "G2M.t0.2", "G3M.t0.3")]) / 3
g_gil_reads <- sum(norm_counts[gil_rows, c("G1.t0.1", "G2.t0.2", "G3.t0.3")]) / 3

gil_fc <- gm_gil_reads / g_gil_reads
gil_fc
norm_counts[gil_rows, c("G1M.t0.1", "G2M.t0.2", "G3M.t0.3", "G1.t0.1", "G2.t0.2", "G3.t0.3")]
