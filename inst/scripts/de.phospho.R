#!/usr/bin/env Rscript
# ==========================
# Phosphopeptide-level differential abundance — CLI entry point
# ==========================
# This is a thin wrapper. All analysis logic lives in R/pipeline.R.
# Load the package, parse arguments, and call run_phospho_pipeline().
#
# During development (before installing):
#   devtools::load_all("/global/projects/proteomics_core/pipelines/proteomics_core_pipelines")
#
# After installing:
#   library(DAPRmd)
#
# Usage:
#   Rscript inst/scripts/de.phospho.R --runid <id> --counts <path> --samplesheet <path> --outdir <dir> --annotation human

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
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Required. Path to merged counts.tsv file"),
  make_option(c("-s", "--samplesheet"), type = "character", default = NULL,
              help = "Required. Path to samplesheet.csv file"),
  make_option(c("-o", "--outdir"), type = "character", default = "./output",
              help = "Output directory [default= %default]"),
  make_option(c("-r", "--runid"), type = "character", default = NULL,
              help = "Required. Unique run ID"),
  make_option(c("-a", "--annotation"), type = "character", default = "mouse",
              help = "Genome for annotation: 'mouse' or 'human' [default= %default]"),
  make_option(c("-i", "--imputation"), type = "character", default = "DEP-MinProb",
              help = paste(
                "Imputation method.",
                "DEP methods (prefix with 'DEP-'): 'DEP-MinProb', 'DEP-knn', 'DEP-bpca', 'DEP-QRILC', 'DEP-man'.",
                "Custom methods: any name matching an impute_<name>() function in R/custom_imputation.R.",
                "'none' skips imputation.",
                "[default= %default]"
              )),
  make_option(c("-q", "--imputation-q"), type = "double", default = 0.01,
              help = "q parameter for MinProb/QRILC: quantile cutoff [default= %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed for reproducibility [default= %default]"),
  make_option(c("--heatmap-top-n"), type = "integer", default = 1000,
              help = "Number of top molecules by CV to show in global heatmap [default= %default]"),
  make_option(c("--heatmap-norm"), type = "character", default = "zscore",
              help = "Heatmap normalization: 'zscore' or 'none' [default= %default]"),
  make_option(c("--gsea-ont"), type = "character", default = "BP",
              help = "Gene ontology category for GSEA: 'BP', 'MF', 'CC', or 'ALL' [default= %default]"),
  make_option(c("--skip-gsea"), action = "store_true", default = FALSE,
              help = "Skip GSEA analysis [default= %default]"),
  make_option(c("--skip-anova"), action = "store_true", default = FALSE,
              help = "Skip one-way ANOVA [default= %default]"),
  make_option(c("--group-color1"), type = "character", default = "#D55E00",
              help = "Color for up-regulated group [default= %default]"),
  make_option(c("--group-color2"), type = "character", default = "#0072B2",
              help = "Color for down-regulated group [default= %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$runid) || is.null(opt$counts) || is.null(opt$samplesheet)) {
  print_help(opt_parser)
  stop("--runid, --counts, and --samplesheet are required.")
}

run_phospho_pipeline(
  run_id            = opt$runid,
  counts_file       = opt$counts,
  samplesheet_file  = opt$samplesheet,
  out_dir           = opt$outdir,
  genome            = opt$annotation,
  imputation_method = opt$imputation,
  imputation_q      = opt$`imputation-q`,
  imputation_seed   = opt$seed,
  heatmap_top_n     = opt$`heatmap-top-n`,
  heatmap_norm      = opt$`heatmap-norm`,
  gsea_ont          = opt$`gsea-ont`,
  skip_gsea         = opt$`skip-gsea`,
  skip_anova        = opt$`skip-anova`,
  group_color1      = opt$`group-color1`,
  group_color2      = opt$`group-color2`
)
