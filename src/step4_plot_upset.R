#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
matrix_path <- NULL
out_pdf <- "consensus.upset.R.pdf"
top_n <- 0
i <- 1
while (i <= length(args)) {
  key <- args[[i]]
  if (key == "--matrix" && i + 1 <= length(args)) {
    matrix_path <- args[[i + 1]]
    i <- i + 2
  } else if (key == "--out-pdf" && i + 1 <= length(args)) {
    out_pdf <- args[[i + 1]]
    i <- i + 2
  } else if (key == "--top-n" && i + 1 <= length(args)) {
    top_n <- as.integer(args[[i + 1]])
    i <- i + 2
  } else if (key == "-h" || key == "--help") {
    cat("Usage: Rscript src/step4_plot_upset.R --matrix <support.tsv> [--out-pdf <file.pdf>] [--top-n <N>]\n")
    quit(status = 0)
  } else {
    stop(paste("Unknown or incomplete argument:", key))
  }
}

if (is.null(matrix_path) || matrix_path == "") {
  stop("--matrix is required.")
}

df <- read.delim(matrix_path, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("chr", "start", "end", "support_count")
if (!all(required_cols %in% colnames(df))) {
  stop("Invalid matrix header. Expected columns: chr, start, end, support_count, <sample IDs...>")
}

sample_ids <- colnames(df)[5:ncol(df)]
if (length(sample_ids) == 0) {
  stop("No sample columns found in matrix.")
}

peak_ids <- paste0(df$chr, ":", df$start, "-", df$end)
peak_list <- lapply(sample_ids, function(sid) {
  keep <- df[[sid]] == 1
  peak_ids[keep]
})
names(peak_list) <- sample_ids

design_mat <- list_to_matrix(peak_list)
comb_mat <- make_comb_mat(design_mat)

if (top_n > 0) {
  ord <- order(comb_size(comb_mat), decreasing = TRUE)
  keep <- ord[seq_len(min(top_n, length(ord)))]
  comb_mat <- comb_mat[, keep]
}

deg <- comb_degree(comb_mat)
palette <- c("#5f75ae", "#ff7f00", "#7ba39d", "#b15928", "#6a3d9a", "#A81E2C")
col_idx <- ((deg - 1) %% length(palette)) + 1
comb_cols <- palette[col_idx]

out_pdf <- normalizePath(out_pdf, winslash = "/", mustWork = FALSE)
out_dir <- dirname(out_pdf)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

pdf(out_pdf, width = 10, height = 7, onefile = TRUE)

p_upset <- UpSet(
  comb_mat,
  set_order = sample_ids,
  comb_col = comb_cols,
  bg_col = "#F0F0FF",
  bg_pt_col = "#CCCCFF",
  top_annotation = upset_top_annotation(
    comb_mat,
    add_numbers = TRUE,
    annotation_name_rot = 90
  ),
  right_annotation = upset_right_annotation(
    comb_mat,
    add_numbers = TRUE,
    show_annotation_name = FALSE,
    axis_param = list(side = "top"),
    gp = gpar(fill = "#08537C", col = "#08537C")
  )
)

draw(p_upset)
dev.off()

cat(sprintf("[OK] Samples: %d\n", length(sample_ids)))
cat(sprintf("[OK] Intersections: %d\n", ncol(comb_mat)))
cat(sprintf("[OK] UpSet PDF written: %s\n", out_pdf))
