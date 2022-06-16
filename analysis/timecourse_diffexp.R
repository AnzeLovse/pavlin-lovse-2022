library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(circlize)
set.seed(123)


source("analysis/helper_functions.R")

RAW_DATA <- file.path("data", "raw")
BASE <- file.path("results", "timecourse_experiment")
dir.create(BASE, showWarnings = FALSE)

counts <- read_counts(path = file.path(RAW_DATA, "timecourse", "merged_counts.tsv"))
coldata <- read_metadata(
  samples = file.path(RAW_DATA, "timecourse", "samples.tsv"),
  units = file.path(RAW_DATA, "timecourse", "units.tsv")
)

counts$G_M.t0.1 <- counts$G.t0.1
counts$G_M.t0.2 <- counts$G.t0.2

zero_rows <- data.frame(rbind(
  c("G_M", "mitomycin", 1, 0),
  c("G_M", "mitomycin", 2, 0)
))

colnames(zero_rows) <- colnames(coldata)
rownames(zero_rows) <- c("G_M.t0.1", "G_M.t0.2")
coldata <- rbind(coldata, zero_rows)

# Filter genes by requiring more than 5 reads in at least two samples.
filtered <- counts[rowSums(counts > 5) >= 2, ]

dds <- DESeqDataSetFromMatrix(
  countData = filtered,
  colData = coldata,
  design = ~ condition + time + condition:time
)
dds <- DESeq(dds, test = "LRT", reduced = ~ condition + time)
resultsNames(dds)

# Plot sample distances
vsd <- vst(dds, blind = TRUE) # variance stabilising transformation
sampleDists <- dist(t(assay(vsd)), method = "euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap_annotation <- HeatmapAnnotation(
  treatment = dds$condition,
  time = dds$time
)

Heatmap(
  sampleDistMatrix,
  col = colors,
  clustering_distance_rows = sampleDists,
  clustering_distance_columns = sampleDists,
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  row_dend_width = unit(30, "mm"),
  top_annotation = heatmap_annotation,
)
ggsave(file.path(BASE, "distance_matrix.pdf"))

# Plot PCA
pca_data <- plotPCA(vsd, intgroup = c("condition", "time"), ntop = nrow(dds), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color = time, shape = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
ggsave(file.path(BASE, "pca.pdf"))

# Extract the LRT test results.
result <- results(dds, alpha = 0.05)
result <- result[order(result$padj), ]

save_deg(
  result = result,
  out_file = file.path(BASE, "timecourse_LRT.tsv"),
  out_annotated = file.path(BASE, "timecourse_LRT_annotated.tsv"),
  fdr_threshold = 0.05
)

comparisons <- list(
  mmc_t5 = "conditionmitomycin.time5",
  mmc_t10 = "conditionmitomycin.time10",
  mmc_t20 = "conditionmitomycin.time20",
  mmc_t30 = "conditionmitomycin.time30"
)

result_list <- list()
shrunken_list <- list()

for (name in names(comparisons)) {
  result <- results(dds, name = comparisons[[name]], test = "Wald", alpha = 0.05)

  shrunken <- lfcShrink(dds, coef = comparisons[[name]], type = "apeglm", res = result)

  result <- result[order(result$padj), ]
  save_deg(
    result = result,
    out_file = file.path(BASE, paste0(name, ".tsv")),
    out_annotated = file.path(BASE, paste0(name, "_annotated.tsv")),
    fdr_threshold = 0.05
  )

  shrunken <- shrunken[order(shrunken$padj), ]
  save_deg(
    result = shrunken,
    out_file = file.path(BASE, paste0(name, "_shrunken.tsv")),
    out_annotated = file.path(BASE, paste0(name, "_shrunken_annotated.tsv")),
    fdr_threshold = 0.05
  )

  result_list[[name]] <- result
  shrunken_list[[name]] <- shrunken
}

# Extract the matrix of log2FC (only maximum likelihood estimates MLE which are not shrunken).
betas <- coef(dds)
res_de <- na.omit(result)[na.omit(result)$padj < 0.05, ]
mat <- betas[rownames(res_de), -c(1, 2)]
mat_shrunken <- as.data.frame(do.call(rbind, lapply(shrunken_list, function(x) x$log2FoldChange)))
col_fun <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

bed <- read.csv(
  file = file.path("data", "external", "annotation_with_gil01.bed"),
  sep = "\t",
  col.names = c("chr", "start", "end", "id", "gene", "strand")
)
rownames(bed) <- bed$id
de_bed <- bed[rownames(mat), ]
de_bed$chr
chromosomes <- bed[rownames(mat), c("chr")]
chr_labes <- as.factor(bed[rownames(mat), "chr"])
levels(chr_labes) <- c("GIL01", "Host", "pBtic")
chr_colors <- setNames(brewer.pal(3, "Pastel2"), levels(chr_labes))
row_annotation <- rowAnnotation(
  chromosome = chr_labes,
  col = list(chromosome = chr_colors)
)

pbtic_deg <- rownames(de_bed[de_bed$chr == "NZ_CP051859.1", ])

pdf(file = file.path(BASE, "pbtic_heatmap.pdf"), width = 12, height = 8)
Heatmap(
  mat[pbtic_deg, 5:8],
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  row_dend_reorder = TRUE,
  row_dend_width = unit(30, "mm"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(title = "log2FC")
)
dev.off()

# Calculate the number of norm counts on left module of gil01.
gene_annotation <- read.csv(
  file = file.path("data/external/annotation_with_gil01.bed"),
  sep = "\t",
  col.names = c("chr", "start", "end", "gene", "type", "strand"),
  header = FALSE
)
gil_ids <- gene_annotation[gene_annotation$chr == "AJ536073.2", "gene"]

norm.count <- counts(dds, normalized = TRUE)
g_m_30 <- sum(norm.count[rownames(norm.count) %in% gil_ids[1:9], c("G_M.t30.1", "G_M.t30.2")]) / 2
g_30 <- sum(norm.count[rownames(norm.count) %in% gil_ids[1:9], c("G.t30.1", "G.t30.2")]) / 2

gil_fc <- g_m_30 / g_30
gil_fc

g_m_20 <- sum(norm.count[rownames(norm.count) %in% gil_ids[1:9], c("G_M.t20.1", "G_M.t20.2")]) / 2
g_20 <- sum(norm.count[rownames(norm.count) %in% gil_ids[1:9], c("G.t20.1", "G.t20.2")]) / 2

gil_fc <- g_m_20 / g_20
gil_fc
