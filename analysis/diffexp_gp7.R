library(DESeq2)
library(ggplot2)
library(RUVSeq)
library(RColorBrewer)
library(EnhancedVolcano)

source("analysis/helper_functions.R")

RAW_DATA <- file.path("data", "raw")
BASE <- file.path("results", "gp7_experiment")
dir.create(BASE, showWarnings = FALSE)

counts <- read_counts(path = file.path(RAW_DATA, "gp7_experiment", "merged_counts.tsv"))
coldata <- read_metadata(
  samples = file.path(RAW_DATA, "gp7_experiment", "samples.tsv"),
  units = file.path(RAW_DATA, "gp7_experiment", "units.tsv")
)

# Remove the two GIL01 samples.
counts <- subset(counts, select = -c(GIL01.t0.1, GIL01_M.t0.1))
coldata <- coldata[!(coldata$sample %in% c("GIL01", "GIL01_M")), ]
coldata$condition <- droplevels(coldata$condition)
counts <- counts[, rownames(coldata)]

# Filter genes by requiring more than 5 reads in at least two samples.
filtered <- counts[rowSums(counts > 5) >= 2, ]

# Plot RLE
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = coldata)
colors <- brewer.pal(6, "Set2")
plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = colors[coldata$condition])
plotRLE(
  betweenLaneNormalization(set, which = "upper"),
  outline = FALSE,
  ylim = c(-4, 4),
  col = colors[coldata$condition]
)

# Get genes that show low differential expression.
pre_dds <- DESeqDataSetFromMatrix(countData = filtered, colData = coldata, design = ~ 0 + condition)
pre_dds <- DESeq(pre_dds)
result_1 <- results(pre_dds, contrast = c("condition", "pDG7_M", "pDG_M"))
result_1 <- result_1[order(result_1$padj), ]

result_2 <- results(pre_dds, contrast = c("condition", "pDG7", "pDG"))
result_2 <- result_2[order(result_2$padj), ]

not_de_1 <- tail(result_1, n = 500)
not_de_2 <- tail(result_2, n = 500)
genes <- intersect(rownames(not_de_1), rownames(not_de_2))

groups <- makeGroups(coldata$condition)
set_norm <- RUVs(x = set, cIdx = genes, k = 2, scIdx = groups)
plotRLE(set_norm, outline = FALSE, ylim = c(-4, 4), col = colors[coldata$condition])

# Run differential expression analysis on the corrected data.

comparisons <- list(
  pDG7_vs_pDG = c("pDG7", "pDG"),
  pDG7_M_vs_pDG_M = c("pDG7_M", "pDG_M"),
  pDG_M_vs_pDG = c("pDG7_M", "pDG_M")
)

labels <- list(
  pDG7_vs_pDG = c(
    "HIS92_RS25225", "HIS92_RS25415", "HIS92_RS24985", # DNA damage
    "HIS92_RS19530", "HIS92_RS22125" # host genes
  ),
  pDG7_M_vs_pDG_M = labels <- c(
    "HIS92_RS25225", "HIS92_RS25415", "HIS92_RS24985", "HIS92_RS15270", "HIS92_RS04600",
    "HIS92_RS06675", # DNA damage
    "HIS92_RS11450", "HIS92_RS17585" # host genes
  ),
  pDG_M_vs_pDG = c(
    "HIS92_RS04150", "HIS92_RS04600", "HIS92_RS06675", "HIS92_RS11675",
    "HIS92_RS11680", "HIS92_RS03705", "HIS92_RS04025", "HIS92_RS15270"
  )
)
result_list <- list()

for (name in names(comparisons)) {
  pData(set_norm)$condition <- relevel(pData(set_norm)$condition, ref = comparisons[[name]][2])
  dds <- DESeqDataSetFromMatrix(
    countData = counts(set_norm),
    colData = pData(set_norm),
    design = ~ W_1 + W_2 + condition
  )

  dds <- DESeq(dds)

  result <- results(dds, contrast = c("condition", comparisons[[name]]), alpha = 0.01)
  condition <- paste(comparisons[[name]], collapse = "_vs_")
  coef <- paste("condition", condition, sep = "_")
  result <- lfcShrink(dds, coef = coef, res = result)
  result <- result[order(result$padj), ]
  
  save_deg(
    result = result,
    out_file = file.path(BASE, paste0(name, ".tsv")),
    out_annotated = file.path(BASE, paste0(name, "_annotated.tsv")),
    fdr_threshold = 0.01
  )

  result_list[[name]] <- result
}

lapply(result_list, function(x) {
  max(abs(x$log2FoldChange))
})
xlim <- c(-4.2, 4.2)
gene_metadata <- get_gene_annot()

for (name in names(result_list)) {
  volcano_title <- stringr::str_replace(string = name, pattern = "_", replacement = " ")
  plot_volcano(
    result = result_list[[name]],
    gene_metadata = gene_metadata,
    outfile = file.path(BASE, paste0(name, "_volcano.pdf")),
    xlim = xlim,
    title = volcano_title,
    fdr_threshold = 0.01,
    labels = labels[[name]]
  )
}
