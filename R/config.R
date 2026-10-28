# ==========================
# Run configuration + comparison definitions
# ==========================
# resolve_run_config()  — merge pipeline defaults < YAML config < CLI flags
# common_cli_options()  — optparse options shared by the inst/scripts/ wrappers
# parse_comparisons()   — turn the `comparisons` setting into contrast definitions
# write_run_config()    — save the fully resolved configuration to the output dir

# Pipeline name -> pipeline function name
PIPELINE_FUNCTIONS <- c(
  proteome = "run_proteome_da_pipeline",
  ptm      = "run_ptm_pipeline",
  ptm_ida  = "run_ptm_ida_pipeline"
)

# CLI flag (optparse dest) -> config key. `--counts` is pipeline-specific and
# handled in cli_to_config().
CLI_TO_CONFIG <- c(
  runid            = "run_id",
  samplesheet      = "samplesheet_file",
  outdir           = "out_dir",
  annotation       = "genome",
  imputation       = "imputation_method",
  `imputation-q`   = "imputation_q",
  seed             = "imputation_seed",
  `heatmap-top-n`  = "heatmap_top_n",
  `heatmap-norm`   = "heatmap_norm",
  `gsea-ont`       = "gsea_ont",
  `skip-gsea`      = "skip_gsea",
  `skip-anova`     = "skip_anova",
  `group-color1`   = "group_color1",
  `group-color2`   = "group_color2",
  `ptm-marker`     = "ptm_marker",
  `protein-counts` = "protein_matrix_file"
)

#' @title Optparse Options Shared by the Pipeline CLI Wrappers
#' @details Every option defaults to \code{NULL} so that \code{\link{resolve_run_config}}
#'   can tell an option given on the command line apart from one left unset; unset
#'   options fall back to the YAML config and then to the pipeline function defaults.
#' @param pipeline Character. One of \code{"proteome"}, \code{"ptm"}, \code{"ptm_ida"}.
#'   Adds the pipeline-specific options (\code{--ptm-marker}, \code{--protein-counts}).
#' @return List of \code{optparse::make_option} objects.
common_cli_options <- function(pipeline = c("proteome", "ptm", "ptm_ida")) {
  pipeline <- match.arg(pipeline)
  counts_help <- if (pipeline == "proteome") {
    "Path to the protein intensity matrix (config key: protein_matrix_file)"
  } else {
    "Path to the PTM peptide intensity matrix (config key: ptm_matrix_file)"
  }

  opts <- list(
    make_option(c("--config"), type = "character", default = NULL,
                help = "YAML run config. Keys are the pipeline argument names; CLI flags override it"),
    make_option(c("-c", "--counts"), type = "character", default = NULL, help = counts_help),
    make_option(c("-s", "--samplesheet"), type = "character", default = NULL,
                help = "Path to samplesheet.csv (SampleID, GroupID)"),
    make_option(c("-o", "--outdir"), type = "character", default = NULL,
                help = "Output directory [default: ./output]"),
    make_option(c("-r", "--runid"), type = "character", default = NULL,
                help = "Unique run ID"),
    make_option(c("-a", "--annotation"), type = "character", default = NULL,
                help = "Genome for annotation: 'mouse' or 'human' [default: mouse]"),
    make_option(c("-i", "--imputation"), type = "character", default = NULL,
                help = paste(
                  "Imputation method.",
                  "DEP methods (prefix with 'DEP-'): 'DEP-MinProb', 'DEP-knn', 'DEP-bpca', 'DEP-QRILC', 'DEP-man'.",
                  "Custom methods: any name matching an impute_<name>() function in R/custom_imputation.R.",
                  "'none' skips imputation.",
                  "[default: DEP-MinProb]"
                )),
    make_option(c("-q", "--imputation-q"), type = "double", default = NULL,
                help = "q parameter for MinProb/QRILC: quantile cutoff [default: 0.01]"),
    make_option(c("--seed"), type = "integer", default = NULL,
                help = "Random seed for reproducibility [default: 42]"),
    make_option(c("--heatmap-top-n"), type = "integer", default = NULL,
                help = "Number of top molecules by CV to show in global heatmap [default: 1000]"),
    make_option(c("--heatmap-norm"), type = "character", default = NULL,
                help = "Heatmap normalization: 'zscore' or 'none' [default: zscore]"),
    make_option(c("--gsea-ont"), type = "character", default = NULL,
                help = "Gene ontology category for GSEA: 'BP', 'MF', 'CC', or 'ALL' [default: BP]"),
    make_option(c("--skip-gsea"), action = "store_true", default = NULL,
                help = "Skip GSEA analysis"),
    make_option(c("--skip-anova"), action = "store_true", default = NULL,
                help = "Skip one-way ANOVA"),
    make_option(c("--group-color1"), type = "character", default = NULL,
                help = "Color for up-regulated group [default: #D55E00]"),
    make_option(c("--group-color2"), type = "character", default = NULL,
                help = "Color for down-regulated group [default: #0072B2]")
  )
  if (pipeline %in% c("ptm", "ptm_ida")) {
    opts <- c(opts, list(
      make_option(c("--ptm-marker"), type = "character", default = NULL,
                  help = "Regex matching PTM peptides in EG.PrecursorId [default: Phospho \\(STY\\)]")
    ))
  }
  if (pipeline == "ptm_ida") {
    opts <- c(opts, list(
      make_option(c("--protein-counts"), type = "character", default = NULL,
                  help = "Path to the imputed protein intensity matrix (config key: protein_matrix_file)")
    ))
  }
  opts
}

#' @title Convert Parsed CLI Options to Config Keys
#' @param opt Named list returned by \code{optparse::parse_args}.
#' @param pipeline Character. One of \code{"proteome"}, \code{"ptm"}, \code{"ptm_ida"}.
#' @return Named list of the options that were given, keyed by pipeline argument name.
cli_to_config <- function(opt, pipeline) {
  key_map <- c(CLI_TO_CONFIG,
               counts = if (pipeline == "proteome") "protein_matrix_file" else "ptm_matrix_file")
  given <- opt[intersect(names(opt), names(key_map))]
  given <- Filter(Negate(is.null), given)
  setNames(given, key_map[names(given)])
}

#' @title Resolve the Configuration for a Pipeline Run
#' @details Merges, in increasing priority: the pipeline function defaults, the YAML
#'   config file, and the CLI options. Unknown keys in the YAML file stop the run so
#'   typos are not silently ignored. Relative paths are resolved against the working
#'   directory, exactly as CLI paths are.
#' @param pipeline Character. One of \code{"proteome"}, \code{"ptm"}, \code{"ptm_ida"}.
#' @param config_file Character or \code{NULL}. Path to a YAML config file.
#' @param cli_args Named list of CLI values keyed by pipeline argument name (see
#'   \code{\link{cli_to_config}}). \code{NULL} entries are ignored.
#' @return Named list of arguments ready for \code{do.call} on the pipeline function.
resolve_run_config <- function(pipeline, config_file = NULL, cli_args = list()) {
  if (!pipeline %in% names(PIPELINE_FUNCTIONS))
    stop(sprintf("Unknown pipeline '%s'. Use one of: %s",
                 pipeline, paste(names(PIPELINE_FUNCTIONS), collapse = ", ")))
  fmls <- formals(get(PIPELINE_FUNCTIONS[[pipeline]]))
  required <- names(fmls)[vapply(fmls, function(x) is.name(x) && !nzchar(as.character(x)), logical(1))]

  cfg <- list()
  if (!is.null(config_file)) {
    if (!file.exists(config_file)) stop(sprintf("Config file not found: %s", config_file))
    cfg <- yaml::read_yaml(config_file)
    if (is.null(cfg)) cfg <- list()
    unknown <- setdiff(names(cfg), names(fmls))
    if (length(unknown) > 0)
      stop(sprintf("Unknown key(s) in %s: %s. Valid keys: %s",
                   config_file, paste(unknown, collapse = ", "), paste(names(fmls), collapse = ", ")))
  }

  cli_args <- Filter(Negate(is.null), cli_args)
  cfg[names(cli_args)] <- cli_args

  # logical defaults (skip_gsea, skip_anova) must stay logical; YAML 'yes'/'no' already are
  for (key in intersect(names(cfg), names(fmls))) {
    if (is.logical(fmls[[key]]) && !(is.logical(cfg[[key]]) && length(cfg[[key]]) == 1))
      stop(sprintf("Config key '%s' must be true or false", key))
  }

  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0)
    stop(sprintf("Missing required setting(s): %s. Provide them in the YAML config or as CLI flags.",
                 paste(missing, collapse = ", ")))
  cfg
}

#' @title Find Legacy 0/1 Comparison Columns in a Samplesheet
#' @param samplesheet Data frame read from the samplesheet CSV.
#' @param exclude Character vector of column names that are never comparison columns.
#' @return Character vector of column names holding only 0, 1 and \code{NA}.
legacy_comparison_columns <- function(samplesheet, exclude = c("SampleID", "GroupID")) {
  candidates <- setdiff(colnames(samplesheet), exclude)
  Filter(function(col) {
    vals <- samplesheet[[col]]
    vals <- vals[!is.na(vals) & vals != ""]
    length(vals) > 0 && all(as.character(vals) %in% c("0", "1"))
  }, candidates)
}

#' @title Build One Comparison Definition
#' @param exp_groups Character vector of experimental (numerator) groups.
#' @param ctrl_groups Character vector of control (denominator) groups.
#' @param name Character or \code{NULL}. Defaults to \code{"<exp>_vs_<ctrl>"}.
#' @return Named list with \code{name}, \code{contrast} (a limma contrast string),
#'   \code{exp}/\code{ctrl} (display labels) and \code{exp_groups}/\code{ctrl_groups}
#'   (groups whose samples are shown in per-comparison plots).
make_comparison <- function(exp_groups, ctrl_groups, name = NULL) {
  group_mean <- function(groups) {
    if (length(groups) == 1) groups
    else sprintf("(%s)/%d", paste(groups, collapse = " + "), length(groups))
  }
  exp_label  <- paste(exp_groups, collapse = "+")
  ctrl_label <- paste(ctrl_groups, collapse = "+")
  list(
    name        = name %||% paste0(exp_label, "_vs_", ctrl_label),
    contrast    = paste(group_mean(exp_groups), "-", group_mean(ctrl_groups)),
    exp         = exp_label,
    ctrl        = ctrl_label,
    exp_groups  = exp_groups,
    ctrl_groups = ctrl_groups
  )
}

#' @title Define the Comparisons for a Run
#' @details Supported values of \code{comparisons}:
#'   \describe{
#'     \item{\code{NULL}}{Use legacy 0/1 comparison columns if the samplesheet has any
#'       (with a deprecation warning), otherwise behave like \code{"all"}.}
#'     \item{\code{"all"}}{Every pair of groups. The group appearing earlier in the
#'       samplesheet is the control, unless \code{reference} is set, in which case the
#'       reference is always the control.}
#'     \item{\code{"vs_reference"}}{Every other group against \code{reference}.}
#'     \item{List of \code{c(exp, ctrl)} pairs}{Exactly these comparisons. The first
#'       group is the numerator, so positive logFC means higher in \code{exp}.}
#'   }
#' @param samplesheet Data frame with \code{SampleID} and \code{GroupID} columns.
#' @param comparisons See Details.
#' @param reference Character or \code{NULL}. Reference (control) group.
#' @return List of comparison definitions as built by \code{\link{make_comparison}}.
parse_comparisons <- function(samplesheet, comparisons = NULL, reference = NULL) {
  groups <- unique(as.character(samplesheet$GroupID))

  bad_names <- groups[make.names(groups) != groups]
  if (length(bad_names) > 0)
    stop(sprintf(paste("GroupID values must be valid R names (letters, digits, '.' or '_',",
                       "not starting with a digit) so limma can build contrasts: %s"),
                 paste(bad_names, collapse = ", ")))

  group_sizes <- table(factor(samplesheet$GroupID, levels = groups))
  small <- names(group_sizes)[group_sizes < 2]
  if (length(small) > 0)
    flog.warn("Group(s) with fewer than 2 samples: %s", paste(small, collapse = ", "))

  if (!is.null(reference) && !reference %in% groups)
    stop(sprintf("reference '%s' is not a GroupID. Groups: %s", reference, paste(groups, collapse = ", ")))

  if (is.null(comparisons)) {
    legacy_cols <- legacy_comparison_columns(samplesheet)
    if (length(legacy_cols) > 0) {
      flog.warn(paste("Comparisons defined by 0/1 samplesheet columns (%s) are deprecated.",
                      "Use `comparisons` in the YAML run config instead."),
                paste(legacy_cols, collapse = ", "))
      result <- list()
      for (col in legacy_cols) {
        vals <- samplesheet[[col]]
        exp_groups  <- unique(as.character(samplesheet$GroupID[vals %in% 1]))
        ctrl_groups <- unique(as.character(samplesheet$GroupID[vals %in% 0]))
        if (length(exp_groups) > 0 && length(ctrl_groups) > 0)
          result[[length(result) + 1]] <- make_comparison(exp_groups, ctrl_groups, name = col)
      }
      return(finalize_comparisons(result))
    }
    comparisons <- "all"
  }

  if (is.character(comparisons) && length(comparisons) == 1 && comparisons == "all") {
    result <- list()
    for (j in seq_along(groups)) {
      for (i in seq_len(j - 1)) {
        pair <- c(groups[j], groups[i])
        if (!is.null(reference) && pair[1] == reference) pair <- rev(pair)
        result[[length(result) + 1]] <- make_comparison(pair[1], pair[2])
      }
    }
    return(finalize_comparisons(result))
  }

  if (is.character(comparisons) && length(comparisons) == 1 && comparisons == "vs_reference") {
    if (is.null(reference)) stop("comparisons: vs_reference requires `reference` to be set")
    result <- lapply(setdiff(groups, reference), function(g) make_comparison(g, reference))
    return(finalize_comparisons(result))
  }

  if (is.character(comparisons) && length(comparisons) == 2) comparisons <- list(comparisons)
  if (!is.list(comparisons))
    stop("`comparisons` must be 'all', 'vs_reference', or a list of [exp, ctrl] pairs")

  result <- lapply(comparisons, function(pair) {
    pair <- unlist(pair)
    if (!is.character(pair) || length(pair) != 2)
      stop(sprintf("Each comparison must be a pair [exp, ctrl]; got: %s",
                   paste(as.character(pair), collapse = ", ")))
    unknown <- setdiff(pair, groups)
    if (length(unknown) > 0)
      stop(sprintf("Comparison [%s] uses unknown group(s): %s. Groups: %s",
                   paste(pair, collapse = ", "), paste(unknown, collapse = ", "),
                   paste(groups, collapse = ", ")))
    if (pair[1] == pair[2]) stop(sprintf("Comparison [%s, %s] compares a group with itself", pair[1], pair[2]))
    make_comparison(pair[1], pair[2])
  })
  finalize_comparisons(result)
}

#' @title Check and Log a Set of Comparison Definitions
#' @param comparisons List of comparison definitions.
#' @return \code{comparisons}, unchanged.
finalize_comparisons <- function(comparisons) {
  if (length(comparisons) == 0) stop("No comparisons defined. Check `comparisons` and the samplesheet GroupID column.")
  comparison_names <- vapply(comparisons, `[[`, character(1), "name")
  dups <- unique(comparison_names[duplicated(comparison_names)])
  if (length(dups) > 0) stop(sprintf("Duplicate comparison(s): %s", paste(dups, collapse = ", ")))
  flog.info("%d comparison(s) defined:\n%s", length(comparisons),
            paste0("  ", comparison_names, ": ", vapply(comparisons, `[[`, character(1), "contrast"),
                   collapse = "\n"))
  comparisons
}

#' @title Write the Resolved Run Configuration
#' @details Writes \code{run_config.yaml} to \code{out_dir} with every pipeline argument
#'   and the comparisons expanded to explicit \code{[exp, ctrl]} pairs, so the run can be
#'   reproduced with \code{--config <out_dir>/run_config.yaml}.
#' @param run_config Named list of the pipeline arguments.
#' @param comparisons List of comparison definitions from \code{\link{parse_comparisons}}.
#' @param out_dir Character. Output directory.
#' @return Invisibly, the path to the written file.
write_run_config <- function(run_config, comparisons, out_dir) {
  run_config$comparisons <- lapply(comparisons, function(cmp) {
    if (length(cmp$exp_groups) == 1 && length(cmp$ctrl_groups) == 1 &&
        cmp$name == paste0(cmp$exp, "_vs_", cmp$ctrl)) {
      c(cmp$exp_groups, cmp$ctrl_groups)
    } else {
      list(name = cmp$name, contrast = cmp$contrast)
    }
  })
  path <- file.path(out_dir, "run_config.yaml")
  writeLines(c(
    sprintf("# Resolved DAPRmd run configuration (%s)", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "# Re-run with: --config run_config.yaml",
    yaml::as.yaml(run_config)
  ), path)
  flog.info("Resolved run configuration written to %s", path)
  invisible(path)
}
