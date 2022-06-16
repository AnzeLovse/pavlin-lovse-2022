library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)

read_counts <- function(path) {
  counts <- read.csv(path, sep = "\t", row.names = "Geneid")
  counts <- subset(counts, select = -c(Chr, Start, End, Strand, Length))

  return(counts)
}

read_metadata <- function(samples, units) {
  samples <- read.csv(file = samples, sep = "\t", stringsAsFactors = TRUE)
  units <- read.csv(file = units, sep = "\t", stringsAsFactors = TRUE)

  coldata <- merge(x = unique(samples), y = units, by = "sample")
  coldata <- subset(x = coldata, select = -c(fq1, fq2))
  coldata[] <- lapply(coldata, as.factor)

  # Sample names are structured as sample_group.time.replicate e.g. pDG7_M.t0.2
  rownames(coldata) <- paste(
    coldata$sample,
    paste0("t", coldata$time),
    coldata$replicate,
    sep = "."
  )

  return(coldata)
}

get_gene_annot <- function() {
  bed <- read.delim(
    file = file.path("data", "external", "annotation_with_gil01.bed"),
    col.names = c("chr", "start", "end", "id", "gene", "strand")
  )
  pbtic_genes <- bed[bed$chr == "NZ_CP051859.1", ]$id
  pbtic_genes_pos <- bed[bed$chr == "NZ_CP051859.1" & bed$strand == "+", ]$id
  pbtic_genes_neg <- bed[bed$chr == "NZ_CP051859.1" & bed$strand == "-", ]$id
  host_genes <- bed[bed$chr == "NZ_CP051858.1", ]$id

  # Cellular response to DNA damage stimulus.
  go <- read.delim(file = file.path("data", "external", "go.tsv"))
  damage_genes <- go[grepl("^GO:0006974$", go$go_id), ]$gene_id

  return(
    list(
      pbtic_genes = pbtic_genes, pbtic_genes_pos = pbtic_genes_pos,
      pbtic_genes_neg = pbtic_genes_neg, host_genes = host_genes,
      damage_genes = damage_genes
    )
  )
}


plot_volcano <- function(result, gene_metadata, outfile, xlim,
                         title, fdr_threshold, labels = c("")) {
  keyvals_colour <- ifelse(
    test = rownames(result) %in% gene_metadata[["damage_genes"]],
    yes = "darkorange",
    no = ifelse(
      test = rownames(result) %in% gene_metadata[["pbtic_genes_pos"]],
      yes = "dodgerblue2",
      no = ifelse(
        test = rownames(result) %in% gene_metadata[["pbtic_genes_neg"]],
        yes = "darkblue",
        no = "black"
      )
    )
  )

  names(keyvals_colour)[keyvals_colour == "darkorange"] <- "Response to DNA damage"
  names(keyvals_colour)[keyvals_colour == "dodgerblue2"] <- "pBtic plasmid genes (+)"
  names(keyvals_colour)[keyvals_colour == "darkblue"] <- "pBtic plasmid genes (-)"
  names(keyvals_colour)[keyvals_colour == "black"] <- "Host genes"

  keyvals_shape <- ifelse(test = result$padj < fdr_threshold, yes = 16, no = 1)
  keyvals_shape[is.na(keyvals_shape)] <- 1
  names(keyvals_shape)[keyvals_shape == "16"] <- paste("p - value <", fdr_threshold)
  names(keyvals_shape)[keyvals_shape == "1"] <- "NS"

  ylim <- c(0, max(-log10(result[["padj"]]), na.rm = TRUE) + 1)

  res_na <- na.omit(result)
  caption_up <- nrow(res_na[res_na$padj < 0.01 & res_na$log2FoldChange > 0, ])
  caption_down <- nrow(res_na[res_na$padj < 0.01 & res_na$log2FoldChange < 0, ])
  caption <- paste(
    "total =", nrow(result), "genes,",
    "up =", caption_up, "genes,",
    "down =", caption_down, "genes"
  )

  other <- subset(rownames(result), !(rownames(result) %in% labels))
  other_labs <- rownames(result)
  other_labs[other_labs %in% other] <- ""

  volcano <- EnhancedVolcano(
    result,
    ylab = bquote(~ -Log[10] ~ italic(Padj)),
    lab = rownames(result),
    title = title,
    subtitle = bquote(italic("Volcano plot (p-adjusted by FC)")),
    x = "log2FoldChange",
    pCutoff = fdr_threshold,
    FCcutoff = 0,
    labSize = 4.0,
    y = "padj",
    colCustom = keyvals_colour,
    shapeCustom = keyvals_shape,
    colAlpha = 0.5,
    selectLab = labels,
    xlim = xlim,
    ylim = ylim,
    drawConnectors = TRUE,
    typeConnectors = "closed",
    arrowheads = FALSE,
    caption = caption,
    max.overlaps = Inf,
    legendPosition = "right"
  )
  ggsave(outfile, volcano, width = 12, height = 8)
}


save_deg <- function(result, out_file, out_annotated, fdr_threshold) {
  out_table <- as.data.frame(result)
  out_table <- cbind(geneid = rownames(out_table), out_table)
  write.table(out_table, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

  gtf <- as.data.frame(
    rtracklayer::import(file.path("data", "external", "annotation_with_gil01.gtf"))
  )
  gene_annotation <- gtf[gtf$type == "gene", c("seqnames", "start", "end", "strand", "gene_id")]
  cds_annotation <- gtf[gtf$type == "CDS", c("gene_id", "gene", "product")]
  transcript_annotation <- gtf[gtf$type == "transcript", c("gene_id", "product")]

  res_na <- na.omit(result)
  res_de <- as.data.frame(res_na[res_na$padj < fdr_threshold, ])

  annotated_out <- merge(
    x = res_de,
    y = gene_annotation,
    by.x = "row.names",
    by.y = "gene_id",
    all.x = TRUE
  )

  annotated_out <- merge(
    x = annotated_out,
    y = cds_annotation,
    by.x = "Row.names",
    by.y = "gene_id",
    all.x = TRUE
  )

  annotated_out <- merge(
    x = annotated_out,
    y = transcript_annotation,
    by.x = "Row.names",
    by.y = "gene_id",
    all.x = TRUE
  )


  # Add tRNA annotation to product x column which has NAs
  no_prod <- is.na(annotated_out$product.x)
  annotated_out[no_prod, ]$product.x <- annotated_out[no_prod, "product.y"]
  annotated_out <- annotated_out[order(annotated_out$padj), ]
  annotated_out <- subset(
    x = annotated_out,
    select = c(
      Row.names, log2FoldChange, lfcSE, padj,
      seqnames, start, end, strand, gene, product.x
    )
  )

  write.table(x = annotated_out, file = out_annotated, sep = "\t", quote = FALSE, row.names = FALSE)
}
