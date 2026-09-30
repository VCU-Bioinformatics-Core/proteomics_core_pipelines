#!/usr/bin/env Rscript
# ==========================
# Single-gene strip plot — CLI entry point
# ==========================
# This is a thin wrapper. Plotting logic lives in R/plot_single_gene.R.
# Plots imputed log2 intensities per group for one gene/protein, with the
# limma statistics from every supplied comparison in a table below.
#
# Usage:
#   Rscript inst/scripts/plot_single_gene.R --gene CD59 \
#     --imputed <outdir>/data/protein_imputed_matrix.csv \
#     --limma <outdir>/data/de_data/compA_limma.csv,<outdir>/data/de_data/compB_limma.csv \
#     --samplesheet <samplesheet.csv> --out CD59_strip_plot.png

library(devtools)
# Resolve the package root from this script's own path, regardless of
# the working directory from which Rscript is invoked.
.script_path <- normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  mustWork = FALSE
)
.pkg_root <- dirname(dirname(dirname(.script_path)))  # inst/scripts/ -> inst/ -> pkg root
load_all(.pkg_root)

library(optparse)

option_list <- list(
  make_option(c("-g", "--gene"), type = "character", default = NULL,
              help = "Required. Gene symbol (e.g. CD59) or UniProt accession (e.g. P13987)"),
  make_option(c("-m", "--imputed"), type = "character", default = NULL,
              help = "Required. Path to imputed matrix CSV (<outdir>/data/protein_imputed_matrix.csv)"),
  make_option(c("-l", "--limma"), type = "character", default = NULL,
              help = "Required. Comma-separated path(s) to <comparison>_limma.csv file(s) (<outdir>/data/de_data/)"),
  make_option(c("-s", "--samplesheet"), type = "character", default = NULL,
              help = "Required. Path to samplesheet.csv file used for the pipeline run"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Output PNG path [default= ./<gene>_strip_plot.png]"),
  make_option(c("--width"), type = "double", default = 7,
              help = "Plot width in inches [default= %default]"),
  make_option(c("--height"), type = "double", default = NULL,
              help = "Plot height in inches [default= 5 + 0.3 per comparison]")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$gene) || is.null(opt$imputed) || is.null(opt$limma) || is.null(opt$samplesheet)) {
  print_help(opt_parser)
  stop("--gene, --imputed, --limma, and --samplesheet are required.")
}

limma_files    <- trimws(strsplit(opt$limma, ",")[[1]])
limma_df       <- read_limma_results(limma_files)
imputed_matrix <- read.csv(opt$imputed, row.names = 1, check.names = TRUE)
samplesheet    <- read.csv(opt$samplesheet)

plots <- plot_single_gene(opt$gene, imputed_matrix, limma_df, samplesheet)

out    <- if (is.null(opt$out)) paste0(opt$gene, "_strip_plot.png") else opt$out
height <- if (is.null(opt$height)) 5 + 0.3 * length(limma_files) else opt$height
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

for (acc in names(plots)) {
  out_file <- if (length(plots) > 1) {
    sub("(\\.[^.]+)$", paste0("_", make.names(acc), "\\1"), out)
  } else out
  save_plot(plots[[acc]], out_file, width = opt$width, height = height)
  flog.info("Saved single-gene plot: %s", out_file)
}
