#!/usr/bin/env Rscript
# ==========================
# Protein-level differential abundance — CLI entry point
# ==========================
# This is a thin wrapper. All analysis logic lives in R/run_proteome_da_pipeline.R.
# Load the package, resolve the run config, and call run_proteome_da_pipeline().
#
# During development (before installing):
#   devtools::load_all("/global/projects/proteomics_core/pipelines/proteomics_core_pipelines")
#
# After installing:
#   library(DAPRmd)
#
# Usage:
#   Rscript inst/scripts/da.proteome.R --config analysis.yaml
#   Rscript inst/scripts/da.proteome.R --runid <id> --counts <path> --samplesheet <path> --outdir <dir> --annotation human
# CLI flags override values from --config.

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

opt_parser <- OptionParser(option_list = common_cli_options("proteome"))
opt        <- parse_args(opt_parser)

cfg <- tryCatch(
  resolve_run_config("proteome", opt$config, cli_to_config(opt, "proteome")),
  error = function(e) {
    print_help(opt_parser)
    stop(conditionMessage(e), call. = FALSE)
  }
)

do.call(run_proteome_da_pipeline, cfg)
