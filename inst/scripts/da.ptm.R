#!/usr/bin/env Rscript
# ==========================
# PTM peptide-level differential abundance — CLI entry point
# ==========================
# This is a thin wrapper. All analysis logic lives in R/run_ptm_pipeline.R.
# Load the package, resolve the run config, and call run_ptm_pipeline().
#
# During development (before installing):
#   devtools::load_all("/global/projects/proteomics_core/pipelines/proteomics_core_pipelines")
#
# After installing:
#   library(DAPRmd)
#
# Usage:
#   Rscript inst/scripts/da.ptm.R --config analysis.yaml
#   Rscript inst/scripts/da.ptm.R --runid <id> --counts <ptm matrix> --samplesheet <path> --outdir <dir> --annotation human
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

opt_parser <- OptionParser(option_list = common_cli_options("ptm"))
opt        <- parse_args(opt_parser)

cfg <- tryCatch(
  resolve_run_config("ptm", opt$config, cli_to_config(opt, "ptm")),
  error = function(e) {
    print_help(opt_parser)
    stop(conditionMessage(e), call. = FALSE)
  }
)

do.call(run_ptm_pipeline, cfg)
