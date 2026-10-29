# ==========================
# Run configuration + comparison definitions
# ==========================
# resolve_run_config()  — merge pipeline defaults < YAML config < CLI flags
# common_cli_options()  — optparse options shared by the inst/scripts/ wrappers
# prepare_samplesheet() — validate samplesheet columns and derive the sample groups
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
                help = "Path to samplesheet.csv (a sample ID column plus variable columns)"),
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

#' @title Validate the Samplesheet and Derive Sample Groups
#' @details The samplesheet has one sample ID column plus any number of variable
#'   columns; the run config says what each variable is for:
#'   \describe{
#'     \item{\code{primary_variables}}{The variables being compared. With one variable,
#'       its values are the groups. With several, e.g. \code{c("tissue", "response")},
#'       each combination is a group named \code{"<tissue>_<response>"} (e.g.
#'       \code{V_R}), giving one limma coefficient per cell of the factorial design.
#'       \code{NULL} uses a column named \code{GroupID} (older two-column samplesheets).}
#'     \item{\code{categorical_covariates}}{Nuisance variables added as categorical
#'       fixed effects (e.g. batch; codes such as 1/2/3 are levels, not numbers).}
#'     \item{\code{numerical_covariates}}{Nuisance variables used as numbers (e.g. age).
#'       They must be numeric with no missing values and not constant; a warning is
#'       logged when the groups explain most of their variation, because the comparisons
#'       then become unreliable.}
#'     \item{\code{blocking_variable}}{One variable whose levels group correlated
#'       samples (e.g. a patient ID).}
#'   }
#'   The rest of the pipeline reads the sample IDs from \code{SampleID} and the
#'   groups from \code{GroupID}, so those columns are (re)written here.
#' @param samplesheet Data frame read from the samplesheet CSV.
#' @param sample_id Character. Column holding the sample IDs (matrix column names).
#' @param primary_variables Character vector of columns, or \code{NULL}.
#' @param categorical_covariates Character vector of columns, or \code{NULL}.
#' @param numerical_covariates Character vector of columns, or \code{NULL}.
#' @param blocking_variable Character or \code{NULL}.
#' @return \code{samplesheet} with \code{SampleID} and \code{GroupID} columns set.
prepare_samplesheet <- function(samplesheet, sample_id = "SampleID", primary_variables = NULL,
                                categorical_covariates = NULL, numerical_covariates = NULL,
                                blocking_variable = NULL) {
  if (!is.character(sample_id) || length(sample_id) != 1)
    stop("`sample_id` must be a single samplesheet column name")
  if (length(blocking_variable) > 1) stop("`blocking_variable` must be a single samplesheet column")
  group_cols <- primary_variables %||% "GroupID"

  missing <- setdiff(c(sample_id, group_cols, categorical_covariates, numerical_covariates, blocking_variable),
                     colnames(samplesheet))
  if (length(missing) > 0) {
    hint <- if (is.null(primary_variables) && "GroupID" %in% missing)
      " Set `primary_variables` to the column(s) that define the groups." else ""
    stop(sprintf("Column(s) not in the samplesheet: %s. Samplesheet columns: %s.%s",
                 paste(missing, collapse = ", "), paste(colnames(samplesheet), collapse = ", "), hint))
  }

  roles <- list(sample_id = sample_id, primary_variables = group_cols,
                categorical_covariates = categorical_covariates, numerical_covariates = numerical_covariates,
                blocking_variable = blocking_variable)
  used <- unlist(roles)
  dups <- unique(used[duplicated(used)])
  if (length(dups) > 0)
    stop(sprintf("Column(s) assigned to more than one of %s: %s",
                 paste(names(roles), collapse = "/"), paste(dups, collapse = ", ")))
  for (col in used) {
    if (any(is.na(samplesheet[[col]]) | samplesheet[[col]] == ""))
      stop(sprintf("Samplesheet column '%s' has empty values", col))
  }
  if (anyDuplicated(samplesheet[[sample_id]]))
    stop(sprintf("Samplesheet column '%s' has duplicated sample IDs", sample_id))
  for (col in numerical_covariates) {
    if (!is.numeric(samplesheet[[col]])) {
      bad <- unique(samplesheet[[col]][is.na(suppressWarnings(as.numeric(as.character(samplesheet[[col]]))))])
      stop(sprintf("numerical_covariates column '%s' must contain only numbers; non-numeric value(s): %s",
                   col, paste(head(bad, 5), collapse = ", ")))
    }
    if (length(unique(samplesheet[[col]])) < 2)
      stop(sprintf("numerical_covariates column '%s' has the same value for every sample", col))
  }

  # internal names: a variable column called SampleID/GroupID would be overwritten
  clash <- intersect(setdiff(used, c(sample_id, group_cols)), c("SampleID", "GroupID"))
  if (sample_id != "SampleID" && "SampleID" %in% group_cols) clash <- c(clash, "SampleID")
  if (length(clash) > 0)
    stop(sprintf("Rename samplesheet column(s) %s: the pipeline uses these names internally",
                 paste(clash, collapse = ", ")))

  samplesheet$SampleID <- as.character(samplesheet[[sample_id]])
  samplesheet$GroupID  <- do.call(paste, c(lapply(samplesheet[group_cols], as.character), sep = "_"))
  group_counts <- table(factor(samplesheet$GroupID, levels = unique(samplesheet$GroupID)))
  flog.info("Groups from %s:\n%s", paste(group_cols, collapse = " x "),
            paste0("  ", names(group_counts), ": ", group_counts, " sample(s)", collapse = "\n"))

  # a numerical covariate that largely tracks the groups (e.g. all responders older)
  # passes the rank check but leaves little information to separate group and covariate
  for (col in numerical_covariates) {
    x <- samplesheet[[col]]
    r2 <- 1 - sum(residuals(lm(x ~ factor(samplesheet$GroupID)))^2) / sum((x - mean(x))^2)
    if (r2 > 0.8)
      flog.warn(paste("The groups explain %.0f%% of the variation in numerical covariate '%s'; group and",
                      "covariate effects are hard to separate, so the comparisons are unreliable"), 100 * r2, col)
  }
  samplesheet
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

#' @title Build a Comparison from a limma Contrast Formula
#' @details Groups with positive coefficients are the experimental side and groups with
#'   negative coefficients the control side; their samples are shown in the
#'   per-comparison plots.
#' @param entry List with \code{name} and \code{contrast}, e.g.
#'   \code{list(name = "R_vs_NR", contrast = "(V_R + L_R)/2 - (V_NR + L_NR)/2")}.
#' @param groups Character vector of all group names.
#' @return Comparison definition (see \code{\link{make_comparison}}).
contrast_comparison <- function(entry, groups) {
  unknown_fields <- setdiff(names(entry), c("name", "contrast"))
  if (length(unknown_fields) > 0)
    stop(sprintf("Unknown field(s) in a contrast entry: %s (use `name` and `contrast`)",
                 paste(unknown_fields, collapse = ", ")))
  contrast <- entry$contrast
  name <- entry$name
  if (!is.character(contrast) || length(contrast) != 1)
    stop("A contrast entry's `contrast` must be a single formula string")
  if (!is.character(name) || length(name) != 1 || !grepl("^[A-Za-z0-9._-]+$", name))
    stop(sprintf("Contrast '%s' needs a `name` made of letters, digits, '.', '_' or '-' (it is used in file names)",
                 contrast))

  coefs <- tryCatch(
    makeContrasts(contrasts = contrast, levels = groups)[, 1],
    error = function(e) stop(sprintf("Invalid contrast '%s' (%s): %s. Groups: %s", name, contrast,
                                     conditionMessage(e), paste(groups, collapse = ", ")), call. = FALSE)
  )
  exp_groups  <- names(coefs)[coefs > 0]
  ctrl_groups <- names(coefs)[coefs < 0]
  if (length(exp_groups) == 0 || length(ctrl_groups) == 0)
    stop(sprintf("Contrast '%s' (%s) needs groups with both positive and negative coefficients", name, contrast))
  if (abs(sum(coefs)) > 1e-8)
    flog.warn("Contrast '%s' coefficients sum to %g, not 0, so it is not a difference between group means",
              name, sum(coefs))
  list(
    name        = name,
    contrast    = contrast,
    exp         = paste(exp_groups, collapse = "+"),
    ctrl        = paste(ctrl_groups, collapse = "+"),
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
#'     \item{List of entries}{Exactly these comparisons. Each entry is either a
#'       \code{c(exp, ctrl)} pair (the first group is the numerator, so positive logFC
#'       means higher in \code{exp}) or a \code{list(name, contrast)} with a limma
#'       contrast formula over the group names, e.g.
#'       \code{"(V_R + L_R)/2 - (V_NR + L_NR)/2"} (see \code{\link{contrast_comparison}}).}
#'   }
#' @param samplesheet Data frame with \code{SampleID} and \code{GroupID} columns (see
#'   \code{\link{prepare_samplesheet}}).
#' @param comparisons See Details.
#' @param reference Character or \code{NULL}. Reference (control) group.
#' @param legacy_columns Logical. Look for legacy 0/1 comparison columns when
#'   \code{comparisons} is \code{NULL}. Only older GroupID samplesheets have them; with
#'   \code{primary_variables} set, 0/1 columns are ordinary variables (e.g. sex).
#' @return List of comparison definitions as built by \code{\link{make_comparison}}.
parse_comparisons <- function(samplesheet, comparisons = NULL, reference = NULL, legacy_columns = TRUE) {
  groups <- unique(as.character(samplesheet$GroupID))

  bad_names <- groups[make.names(groups) != groups]
  if (length(bad_names) > 0)
    stop(sprintf(paste("Group names (from primary_variables) must be valid R names (letters, digits, '.' or '_',",
                       "not starting with a digit) so limma can build contrasts: %s"),
                 paste(bad_names, collapse = ", ")))

  group_sizes <- table(factor(samplesheet$GroupID, levels = groups))
  small <- names(group_sizes)[group_sizes < 2]
  if (length(small) > 0)
    flog.warn("Group(s) with fewer than 2 samples: %s", paste(small, collapse = ", "))

  if (!is.null(reference) && !reference %in% groups)
    stop(sprintf("reference '%s' is not a group. Groups: %s", reference, paste(groups, collapse = ", ")))

  if (is.null(comparisons)) {
    legacy_cols <- if (legacy_columns) legacy_comparison_columns(samplesheet) else character(0)
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
  if (is.list(comparisons) && !is.null(comparisons$contrast)) comparisons <- list(comparisons)
  if (!is.list(comparisons))
    stop("`comparisons` must be 'all', 'vs_reference', or a list of [exp, ctrl] pairs and {name, contrast} entries")

  result <- lapply(comparisons, function(entry) {
    if (is.list(entry) && !is.null(names(entry))) return(contrast_comparison(entry, groups))
    pair <- unlist(entry)
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
  if (length(comparisons) == 0) stop("No comparisons defined. Check `comparisons` and `primary_variables`.")
  comparison_names <- vapply(comparisons, `[[`, character(1), "name")
  dups <- unique(comparison_names[duplicated(comparison_names)])
  if (length(dups) > 0) stop(sprintf("Duplicate comparison(s): %s", paste(dups, collapse = ", ")))
  flog.info("%d comparison(s) defined:\n%s", length(comparisons),
            paste0("  ", comparison_names, ": ", vapply(comparisons, `[[`, character(1), "contrast"),
                   collapse = "\n"))
  comparisons
}

#' @title Stop Unless Every Comparison Is a Simple Two-Group Pair
#' @details Used by pipelines that test each comparison with per-feature t-tests rather
#'   than a limma contrast (PTM-IDA), where pooled groups and weighted contrasts have no
#'   equivalent.
#' @param comparisons List of comparison definitions.
#' @param pipeline_label Character. Pipeline name for the error message.
#' @return \code{TRUE} invisibly.
check_simple_comparisons <- function(comparisons, pipeline_label) {
  is_simple <- vapply(comparisons, function(cmp) {
    length(cmp$exp_groups) == 1 && length(cmp$ctrl_groups) == 1 &&
      cmp$contrast == paste(cmp$exp_groups, "-", cmp$ctrl_groups)
  }, logical(1))
  if (!all(is_simple))
    stop(sprintf("%s supports only [exp, ctrl] pair comparisons. Not supported: %s", pipeline_label,
                 paste(vapply(comparisons[!is_simple], function(cmp) sprintf("%s (%s)", cmp$name, cmp$contrast),
                              character(1)), collapse = "; ")))
  invisible(TRUE)
}

#' @title Write the Resolved Run Configuration
#' @details Writes \code{run_config.yaml} to \code{out_dir} with every pipeline argument
#'   and the comparisons expanded to explicit \code{[exp, ctrl]} pairs (or
#'   \code{name}/\code{contrast} entries), so the run can be reproduced with
#'   \code{--config <out_dir>/run_config.yaml}.
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
