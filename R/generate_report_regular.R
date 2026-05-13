
# source("report_generator.R"); generate_report('analysis.rds', output_dir = getwd()) # test below
generate_report_regular <- function(analysis_results_path, output_dir = "./", report_prefix = "proteomics_analysis", analyst="Joaquin Reyna") {
  
  # Validate inputs
  if (!file.exists(analysis_results_path)) {
    stop("Analysis results file does not exist:", analysis_results_path)
  }

  # Resolve to absolute paths so rmarkdown::render() can find them
  # regardless of working directory changes during rendering
  analysis_results_path <- normalizePath(analysis_results_path, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate unique filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  #output_file <- file.path(output_dir, paste0(report_prefix, "_", timestamp, ".html"))
  output_file <- file.path(output_dir, glue("{report_prefix}.html"))
  
  # Load results
  rds_data <- readRDS(analysis_results_path)
  
  # Extract components from RDS
  results <- rds_data[[1]]
  comparisons <- rds_data[[2]]
  out_dirs <- rds_data[[3]]
  pca_plot <- rds_data[[4]]
  #pca_plotly <- rds_data[[5]]
  #pca_3d <- rds_data[[6]]
  #annotation <- rds_data[[7]]
  
  # calculate the current date
  date <- format(Sys.time(), "%B %d, %Y")
  
  # Create R Markdown template (beginning section)
  rmd_content <- glue('---
title: "Proteomics Differential Abundance Analysis Report"
author: "Bioinformatics Shared Resources at VCU"
date: "{date}"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
    highlight: tango
    df_print: paged
    code_folding: hide
---

```{{r setup, include=FALSE}}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE,
                      dpi = 150, fig.retina = 2)









# Load results
rds_data <- readRDS("{analysis_results_path}")
results <- rds_data[[1]]
comparisons <- rds_data[[2]]
out_dirs <- rds_data[[3]]
pca_plot <- rds_data[[4]]
#pca_plotly <- rds_data[[5]]
intensity_matrix_raw <- rds_data[[5]]
intensity_matrix <- rds_data[[6]]
#annotation <- rds_data[[7]]
sample_info <- rds_data[[8]]
protein_counts <- rds_data[[9]]
analysis_params <- rds_data[[10]]
anova_summary <- rds_data[[11]]
color1 <- if (!is.null(analysis_params$color1)) analysis_params$color1 else "#D55E00"
color2 <- if (!is.null(analysis_params$color2)) analysis_params$color2 else "#0072B2"
```

## Overview

**Analyst: {analyst}**

**Report Guide:** This report details the analysis pipeline, intermediate processing steps, and
final outputs. We begin with **Global Sample Exploration**, which covers protein-level
normalization, imputation, and sample relationships visualized through PCA and global heatmaps.
The **Differential Abundance Results** follow, starting with a master "Query Differential Expression Hits Across All Analyses
" table for a quick
cross-comparison of significant hits. Detailed subsections for each comparison provide volcano
plots, MA plots, heatmaps, and comprehensive protein tables. For deeper dives into the raw data
or downstream Ingenuity Pathway Analysis (IPA), please see the output directories or the
**Utilizing IPA** section. Suggested draft language can be found in **Manuscript-Ready Text**.

**Brief synopsis:** This analysis report contains the results for `r length(comparisons)` comparisons using differential abundance analysis.


## Pipeline

For this analysis we used the following steps:

1. **Imputation**: Missing values were imputed prior to differential abundance testing.
2. **Differential abundance**: Limma was used to identify differentially abundant proteins.
3. **Functional analysis**: Gene Set Enrichment Analysis (GSEA) was performed using clusterProfiler.
4. **Thresholds**: Proteins with adjusted p-value < 0.05 and |FC| >= 1.5 (equivalent to log2(FC) = 0.58) were considered differentially abundant.

As part of this pipeline we produce the following files for your downstream use:

```
.
├── data
│   ├── anova
│   │   └── global_anova.csv
│   ├── de_data
│   │   └── [comparison]_limma.csv
│   └── gsea_data
│       └── [comparison]_go_analysis_.csv
└── figures
    ├── gsea
    │   └── [comparison]_gsea.png
    ├── heatmap
    │   ├── global_heatmap.png
    │   └── [comparison]_heatmap.png
    ├── ma
    │   └── [comparison]_ma.png
    ├── pca
    │   └── global_pca_plot.png
    └── volcano
        └── [comparison]_volcano.png
```

### Run Summary

```{{r run-summary}}
analyses <- data.frame(
  Level = c("Global", "Per-Comparison", "Per-Comparison"),
  Analysis = c("One-Way ANOVA", "Differential Abundance (Limma)", "GSEA"),
  Status = c(
    ifelse(isTRUE(anova_summary$skipped), "\U1F6AB",
           ifelse(anova_summary$n_groups > 2, "\u2705", "N/A")),
    "\u2705",
    ifelse(isTRUE(analysis_params$skip_gsea), "\U1F6AB", "\u2705")
  ),
  stringsAsFactors = FALSE
)
knitr::kable(analyses, col.names = c("Level", "Step", "Status"))
```

### Analysis Parameters

```{{r analysis-parameters}}
if (!is.null(analysis_params)) {{
  ont_labels <- c(
    BP  = "Biological Process (BP)",
    MF  = "Molecular Function (MF)",
    CC  = "Cellular Component (CC)",
    ALL = "All ontologies (ALL)"
  )
  ont_label <- ifelse(
    !is.null(analysis_params$gsea_ont) &&
      analysis_params$gsea_ont %in% names(ont_labels),
    ont_labels[[analysis_params$gsea_ont]],
    analysis_params$gsea_ont
  )
  gsea_performed <- ifelse(isTRUE(analysis_params$skip_gsea), "No", "Yes")
  params_df <- data.frame(
    Parameter = c(
      "Genome / Annotation",
      "GSEA Gene Ontology",
      "GSEA Performed",
      "Global Heatmap Top-N"
    ),
    Value = c(
      tools::toTitleCase(analysis_params$genome),
      ont_label,
      gsea_performed,
      as.character(analysis_params$heatmap_top_n)
    ),
    stringsAsFactors = FALSE
  )
  knitr::kable(params_df)
}}
```

## Global Sample Exploration

### Protein Filtering Summary

The following table tracks how many proteins remain at each stage of the filtering pipeline — from the initial search results through contaminant (cRAP) removal and imputation eligibility checks. Only proteins that pass all filters are carried forward into normalization and differential abundance testing.

```{{r global-summary}}
if (!is.null(protein_counts)) {{
  steps  <- c("Total proteins", "After cRAP removal")
  counts <- c(protein_counts$total, protein_counts$no_crap)
  if (!is.null(protein_counts$not_imputable) && protein_counts$not_imputable > 0) {{
    steps  <- c(steps,  "Not-imputable (discarded)")
    counts <- c(counts, protein_counts$not_imputable)
    steps  <- c(steps,  "Curated proteins (used downstream)")
    counts <- c(counts, protein_counts$no_crap - protein_counts$not_imputable)
  }}
  summary_df <- data.frame(Step = steps, Count = counts, stringsAsFactors = FALSE)
  knitr::kable(summary_df)
}}
```
\n')

rmd_content <- paste0(rmd_content, '
### Protein Levels Across Observed Values Per Sample

These boxplots display the distribution of raw (pre-imputation) observed intensities for each sample. Only measured values are included — missing values are excluded. Comparable median intensities and spread across samples indicate consistent sample loading and data quality, while large deviations may suggest normalization issues or outlier samples.

```{r protein-boxplot, fig.width=max(8, ncol(intensity_matrix_raw) * 0.6), fig.height=5}
# Reshape raw matrix (pre-imputation) to long format, keeping only observed (non-NA) values
bp_long <- as.data.frame(intensity_matrix_raw) %>%
  tibble::rownames_to_column("protein") %>%
  tidyr::pivot_longer(-protein, names_to = "sample", values_to = "intensity") %>%
  dplyr::filter(!is.na(intensity))

# Preserve sample order from the matrix columns
bp_long$sample <- factor(bp_long$sample, levels = colnames(intensity_matrix_raw))

ggplot(bp_long, aes(x = sample, y = intensity)) +
  geom_boxplot(fill = color2, outlier.size = 0.8, outlier.alpha = 0.5) +
  labs(title = "Protein Levels Across Observed Values Per Sample",
       x = "Sample", y = "Intensity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
```

')

imputation_content <- '\n
### Imputation Reporting

This section documents the imputation strategy applied to handle missing values in the dataset. The parameters used are recorded for reproducibility, followed by a series of diagnostic plots that characterize the extent and distribution of imputed values across samples and proteins.

### Parameters

```{r imputation-params}
imp_params <- rds_data[[7]]
imp_method <- if (!is.null(imp_params)) imp_params$method else "unknown"
imp_q      <- if (!is.null(imp_params)) imp_params$q      else NA

param_table <- data.frame(
  Parameter = c("Package", "Method", "q (quantile cutoff)"),
  Value     = c(if (startsWith(imp_method, "DEP-")) "DEP" else "VCU Core", imp_method, as.character(imp_q)),
  stringsAsFactors = FALSE
)
kable(param_table)
```

### Observed vs. Imputed Intensity Values

The histogram below shows the distribution of all log2 intensity values across
all proteins. **Observed** (blue) values are those measured directly, while
**Imputed** (orange) values are those that were missing and filled in. Imputed
values typically appear as a secondary peak shifted towards the lower end of
the distribution, reflecting the left-censored nature of missing-not-at-random
(MNAR) data in proteomics.

```{r imputation-histogram, out.width="100%"}
knitr::include_graphics(file.path(out_dirs$imputation, "global_imputation_histogram.png"))
```

### Observed vs. Imputed Value Counts

The table below summarizes the total number of observed and imputed intensity values across all curated proteins and samples. Observed values were measured directly by the instrument; imputed values replaced missing entries using the selected imputation method. A high imputation percentage warrants caution — it suggests a substantial portion of the data is modeled rather than measured.

```{r imputation-counts}
intensity_raw     <- rds_data[[5]]
intensity_imputed <- rds_data[[6]]
raw_mat     <- as.matrix(intensity_raw)
imputed_mat <- as.matrix(intensity_imputed)
raw_mat  <- raw_mat[rownames(imputed_mat), , drop = FALSE]
obs_mask <- !is.na(raw_mat)
n_obs    <- sum(obs_mask)
n_imp    <- sum(!obs_mask)
n_tot    <- length(obs_mask)
pct_imp  <- round(100 * n_imp / n_tot, 1)
n_curated <- nrow(intensity_imputed)
counts_table <- data.frame(
  Parameter = c("Curated proteins", "Total measurements", "Observed values", "Imputed values", "Imputed (%)"),
  Value     = c(format(n_curated, big.mark = ","),
                format(n_tot,     big.mark = ","),
                format(n_obs,     big.mark = ","),
                format(n_imp,     big.mark = ","),
                paste0(pct_imp, "%")),
  stringsAsFactors = FALSE
)
kable(counts_table)
```

### Number of Missing Values per Sample

This bar chart shows how many missing values were present in each sample before imputation. Samples with a disproportionately high number of missing values may indicate poor sample quality or technical issues and should be reviewed carefully.

```{r imputation-missing-per-sample, out.width="100%"}
knitr::include_graphics(file.path(out_dirs$imputation, "global_imputation_missing_per_sample.png"))
```

### Total Imputed Values per Protein

This chart summarizes how many values were imputed for each protein across all samples. Proteins with a high proportion of imputed values should be interpreted with caution, as their quantification relies heavily on modeled rather than measured intensities.

```{r imputation-total-counts, out.width="100%"}
knitr::include_graphics(file.path(out_dirs$imputation, "global_imputation_total_counts.png"))
```

### Distribution of Imputed Values per Protein

This plot shows the per-protein distribution of imputed value counts. It helps identify whether imputation is concentrated in a small subset of proteins or spread broadly across the dataset, informing confidence in downstream differential abundance results.

```{r imputation-per-protein, out.width="100%"}
knitr::include_graphics(file.path(out_dirs$imputation, "global_imputation_distribution.png"))
```
'

rmd_content <- paste0(rmd_content, imputation_content)

rmd_content <- paste0(rmd_content, glue('\n
### Principal Component Analysis

PCA was performed to visualize the overall patterns of protein levels across samples and
to identify potential batch effects or outliers. **What we expect:** Samples with similar
profiles should cluster together, while dissimilar samples are expected to
separate into distinct clusters.

```{{r pca-plot, fig.width=7, fig.height=6 }}
pca_plot
```

```{{r anova-analysis, results="asis"}}
if (!isTRUE(anova_summary$skipped)) {{
  cat("### One-Way ANOVA\n\n")
  if (anova_summary$n_groups > 2) {{
    anova_path <- file.path(out_dirs$anova, "global_anova.csv")
    if (file.exists(anova_path)) {{
      anova_df <- read.csv(anova_path)
      cat("A one-way ANOVA was performed across all", anova_summary$n_groups, "groups to identify",
          "proteins that vary significantly beyond any single pairwise comparison.",
          "P-values were adjusted using the Benjamini-Hochberg (BH) method.\n\n")
      cat(paste0("**", anova_summary$n_sig, "** out of **", anova_summary$n_total,
                 "** proteins are significantly variable across groups",
                 " (ANOVA FDR < 0.05).\n\n"))
      print(knitr::kable(head(anova_df, 20),
        caption = "Top 20 proteins by one-way ANOVA (sorted by adjusted p-value)",
        row.names = FALSE))
    }}
  }} else {{
    cat("*One-way ANOVA is only reported when more than two groups are present.*\n\n")
  }}
}}
```

### Global Protein Heatmap

```{{r global-heatmap-desc, results="asis"}}
if (isTRUE(analysis_params$heatmap_norm == "zscore")) {{
  cat("The heatmap below shows z-score normalized expression of the top 1000 proteins
by coefficient of variation (CV) across all samples, clustered by similarity.")
}} else {{
  cat("The heatmap below shows intensity expression of the top 1000 proteins
by coefficient of variation (CV) across all samples, clustered by similarity.")
}}
```

```{{r global-heatmap, out.width="100%"}}
global_heatmap_path <- file.path(out_dirs$heatmap, "global_heatmap.png")
if (file.exists(global_heatmap_path)) {{
  knitr::include_graphics(global_heatmap_path)
}} else {{
  cat("Global heatmap not available.")
}}
```

<!--Interactive PCA plots can also be found in the following files:

- 2D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot.html")`
- 3D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot3D.html")`-->

```{{r results-summary}}
# Create a summary table for all comparisons
n_ensembl_mapped <- if (!is.null(protein_counts$ensembl_mapped)) protein_counts$ensembl_mapped else NA_integer_

summary_table <- data.frame(
  Comparison = character(),
  Experimental = character(),
  Control = character(),
  Total_DAPs = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  `# Ensembl Mapped` = integer(),
  GSEA_Performed = character(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

for (i in seq_along(comparisons)) {{
  if (i <= length(results) && !is.null(results[[i]])) {{
    res_df <- results[[i]]$limma

    # Count DAPs (padj < 0.05 & |log2FC| >= 0.58)
    if (!is.null(res_df)) {{
      total_daps <- sum(!is.na(res_df$adj.P.Val) & res_df$adj.P.Val < 0.05 & abs(res_df$logFC) >= 0.58)
      up_daps <- sum(!is.na(res_df$adj.P.Val) & res_df$adj.P.Val < 0.05 & res_df$logFC >= 0.58)
      down_daps <- sum(!is.na(res_df$adj.P.Val) & res_df$adj.P.Val < 0.05 & res_df$logFC <= -0.58)

      gsea_status <- ifelse(!is.null(results[[i]]$gsea), "Yes", "No")

      summary_table <- rbind(summary_table, data.frame(
        Comparison = comparisons[[i]]$name,
        Experimental = comparisons[[i]]$exp,
        Control = comparisons[[i]]$ctrl,
        Total_DAPs = total_daps,
        Upregulated = up_daps,
        Downregulated = down_daps,
        `# Ensembl Mapped` = n_ensembl_mapped,
        GSEA_Performed = gsea_status,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
    }}
  }}
}}
```
  
## Differential Abundance Results
This section contains a high-level table summarizing the differential abundance results,
followed by a subsection for each comparison that visualizes the results using a volcano
plot, heatmap, table of top hits, and GSEA results (shown as both a bubble plot and a table).

### Summary of Comparisions

Briefly, we report the following:

- Total number of comparisons analyzed: `r length(comparisons)`
- Comparisons with the highest number of DAPs: `r summary_table$Comparison[which.max(summary_table$Total_DAPs)]` (`r max(summary_table$Total_DAPs)` DAPs)
- Comparisons with the lowest number of DAPs: `r summary_table$Comparison[which.min(summary_table$Total_DAPs)]` (`r min(summary_table$Total_DAPs)` DAPs)

The table below summarizes differential abundance results from all comparisons.
Specifically, the **Total_DAPs** column reports the number of proteins with an adjusted
p-value < 0.05 and |FC| ≥ 1.5. The **Upregulated** and **Downregulated** columns break
down this count accordingly. The **GSEA_Performed** column indicates whether GSEA was run
for each case.

```{{r display-summary-table}}
# Display the summary table
kable(summary_table, caption = "")
```

### Query Differential Expression Hits Across All Analyses

Use the search boxes below to look up any protein. **Protein Name** is the gene symbol; **UniProt ID** is the Swiss-Prot accession. **Significant Comparisons** lists every comparison in which the protein was significant (adj. p-value < 0.05 and |FC| ≥ 1.5), or "Not Significant" if it did not meet that threshold in any comparison. Use the column search boxes (row below the header) to filter by protein name or UniProt ID — each column search box filters independently.

```{{r query-protein}}
all_rows <- dplyr::bind_rows(lapply(seq_along(comparisons), function(i) {{
  res <- results[[i]]$limma
  if (is.null(res)) return(NULL)
  res$comparison_name <- comparisons[[i]]$name
  res$is_sig <- !is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & abs(res$logFC) >= 0.58
  res[, intersect(c("gene_name", "uniprotswissprot", "comparison_name", "is_sig"), colnames(res))]
}}))
id_cols <- intersect(c("gene_name", "uniprotswissprot"), colnames(all_rows))
if (nrow(all_rows) == 0 || length(id_cols) == 0) {{
  query_df <- data.frame()
}} else {{
  query_df <- all_rows %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
    dplyr::summarise(
      Significant_Comparisons = {{
        sig_comps <- comparison_name[is_sig]
        if (length(sig_comps) == 0) "Not Significant" else paste(sig_comps, collapse = "<br>")
      }},
      .groups = "drop"
    )
  if ("gene_name" %in% colnames(query_df)) {{
    query_df <- query_df %>% dplyr::arrange(Significant_Comparisons == "Not Significant", gene_name)
  }} else {{
    query_df <- query_df %>% dplyr::arrange(Significant_Comparisons == "Not Significant")
  }}
  query_df <- query_df %>%
    dplyr::rename(dplyr::any_of(c("Protein Name" = "gene_name", "UniProt ID" = "uniprotswissprot",
                                  "Significant Comparisons" = "Significant_Comparisons")))
}}
DT::datatable(query_df,
              caption = "Protein significance across all comparisons",
              rownames = FALSE,
              escape = FALSE,
              filter = "top",
              options = list(pageLength = 15, dom = "tip"))
```
\n'))
  
  # Add detailed sections for each comparison
  for (i in seq_along(comparisons)) {
    
    name = comparisons[[i]]$name
    exp = comparisons[[i]]$exp
    ctrl = comparisons[[i]]$ctrl
    
    comparison_section <- glue('\n
### {name}

**Experimental design**

- **Experimental group**: {exp}
- **Control group**: {ctrl}


**Volcano Plot**

- Description: main visualization for differential abundance results. This rendition uses
  an **orange horizontal line** to indicate the significant p-value threshold of \\< 0.05. Therefore,
  every point (protein) above that orange line can be considered statistically significant (note:
  before FDR correction). In addition, the **black vertical lines** indicate 1.5 fold change.
  Proteins highlighted in **orange** are more abundant in **{exp}** relative to **{ctrl}** (positive logFC).
  Proteins highlighted in **blue** are more abundant in **{ctrl}** relative to **{exp}** (negative logFC).
  Proteins in gray do not meet the thresholds for both logFC and p-value.
- Data point: a protein
- X-axis: log2 fold-change ({exp} / {ctrl}) — positive = higher in {exp}
- Y-axis: -log10(p-value)

```{{r volcano-{i}, out.width="90%" }}
# Display volcano plot from file
volcano_path <- file.path(out_dirs$volcano, paste0("{name}_volcano.png"))
if (file.exists(volcano_path)) {{
  knitr::include_graphics(volcano_path)
}} else {{
  cat("Volcano plot not available for this comparison")
}}
```


**MA Plot**

- Description: visualizes the relationship between average expression and fold change.
  Points are colored by significance status. A horizontal dashed line marks y = 0
  (no change), and dotted lines mark ±log2(1.5) fold change thresholds.
- Data point: a protein
- X-axis: average log2 intensity across all samples (AveExpr)
- Y-axis: log2 fold-change ({exp} / {ctrl})

```{{r ma-plot-{i}, out.width="90%"}}
ma_path <- file.path(out_dirs$ma, paste0("{name}_ma.png"))
if (file.exists(ma_path)) {{
  knitr::include_graphics(ma_path)
}} else {{
  cat("MA plot not available for this comparison")
}}
```

**Heatmap**

```{{r heatmap-desc-{i}, results="asis"}}
if (isTRUE(analysis_params$heatmap_norm == "zscore")) {{
  cat("- Description: a heatmap of **z-score normalized** intensity data with application of **hierarchical clustering** by both samples and protein levels\n")
  cat("- X-axis: samples\n")
  cat("- Y-axis: proteins\n")
  cat("- Color-scale: z-score normalized intensity levels\n")
}} else {{
  cat("- Description: a heatmap of **intensity** data with application of **hierarchical clustering** by both samples and protein levels\n")
  cat("- X-axis: samples\n")
  cat("- Y-axis: proteins\n")
  cat("- Color-scale: intensity levels\n")
}}
```


```{{r heatmap-{i}}}
# Display heatmap from file
heatmap_path <- file.path(out_dirs$heatmap, paste0("{name}_heatmap.png"))
if (file.exists(heatmap_path)) {{
  knitr::include_graphics(heatmap_path)
}} else {{
  cat("Heatmap not available for this comparison")
}}
```


#### Significant DAPs by Imputation Category

```{{r top-daps-{i} }}
# Display top DAPs table

dap_flag = 0
if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$limma)) {{
  imp_levels <- c("complete-data", "on-off", "imputation-low", "imputation-medium", "imputation-high", "not-significant", "other")

  combined_sigs <- results[[{i}]]$limma %>%
    filter(!is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= 0.58) %>%
    mutate(logFC = round(logFC, 2),
           pvalue = signif(P.Value, 3),
           padj = signif(adj.P.Val, 3),
           imputation_category = factor(imputation_category, levels = imp_levels)) %>%
    arrange(padj)

  # update the flags
  if (nrow(combined_sigs) > 0) {{
    dap_flag <- 1
  }}

}}
```

This table breaks down all significant proteins (adj. p-value < 0.05 and |FC| ≥ 1.5) by their imputation category, giving an overall picture of how much the significant hits relied on imputed values. Proteins classified as *complete-data* had no missing values, while *imputation-low/medium/high* indicate increasing levels of imputation. *on-off* is a special case where the protein is fully observed in one group but completely missing in the other, representing a presence/absence pattern rather than a continuous abundance difference.

```{{r imputation-category-summary-{i} }}
if (i <= length(results) && !is.null(results[[{i}]]) && !is.null(results[[{i}]]$limma)) {{
  sig_df <- results[[{i}]]$limma %>%
    filter(!is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= 0.58)
  if ("imputation_category" %in% colnames(sig_df) && nrow(sig_df) > 0) {{
    cat_counts <- as.data.frame(table(sig_df$imputation_category))
    colnames(cat_counts) <- c("Category", "Count")
    kable(cat_counts)
  }}
}}
```

#### Table of Differentially Abundant Proteins

All significant proteins (adj. p-value < 0.05 and |FC| ≥ 1.5), sorted by adjusted p-value.

```{{r daps-desc-{i}, results="asis"}}

if (dap_flag > 0){{
  cat("- Description: all significant differentially abundant proteins from the limma moderated t-test\n")
  cat("- **Protein Name**: human-readable protein name\n")
  cat("- **UniProt ID**: UniProt accession identifier\n")
  cat("- **logFC**: log\u2082 fold-change (experimental / control) — positive = higher in experimental, negative = higher in control\n")
  cat("- **P-value**: raw p-value from the limma moderated t-test\n")
  cat("- **Adjusted P-value**: Benjamini-Hochberg FDR-adjusted p-value\n")
  cat("- **Imputation Category**: extent of missing value imputation — *complete-data* (none), *imputation-low* (1 missing), *imputation-medium* (2 missing), *imputation-high* (3+ missing)\n\n")
}}

```

```{{r daps-{i}}}

if (dap_flag > 0){{
  DT::datatable(
    combined_sigs %>%
      dplyr::select(dplyr::any_of(c("gene_name", "uniprotswissprot")),
                    logFC, pvalue, padj, imputation_category) %>%
      dplyr::rename(dplyr::any_of(c("Protein Name"       = "gene_name",
                                    "UniProt ID"         = "uniprotswissprot",
                                    "P-value"            = "pvalue",
                                    "Adjusted P-value"   = "padj",
                                    "Imputation Category" = "imputation_category"))),
    caption = "All Significant Differentially Abundant Proteins",
    filter = "top")

}} else {{
  cat("No differential abundance results available for this comparison\n")
}}
  
```

```{{r gsea-header-{i}, results="asis"}}
if (!isTRUE(analysis_params$skip_gsea)) {{
  cat("**Gene Set Enrichment Analysis**\n\n")
  cat("- Description: Gene Set Enrichment Analysis using Gene Ontology (GO) terms, displayed as a bar plot and table\n")
  cat("- X-axis: Normalized Enrichment Score (NES) — a positive NES indicates that the gene set is activated (upregulated) in the experimental condition, while a negative NES indicates that the gene set is suppressed (downregulated)\n")
  cat("- Y-axis: a given Gene Set\n\n")
}}
```

```{{r gsea-{i}, out.width="100%"}}
# Display GSEA results from file
if (!isTRUE(analysis_params$skip_gsea)) {{
  gsea_path <- file.path(out_dirs$gsea, paste0("{name}_gsea.png"))
  if (file.exists(gsea_path)) {{
    knitr::include_graphics(gsea_path)
  }}
}}
```

```{{r gsea-desc-{i}, results="asis"}}
if (!isTRUE(analysis_params$skip_gsea)) {{
  if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$gsea)) {{
    gsea_df <- as.data.frame(results[[{i}]]$gsea)
    if (nrow(gsea_df) > 0) {{
      cat("- **ID**: Gene Ontology term identifier\n")
      cat("- **Description**: Gene Ontology term name\n")
      cat("- **setSize**: number of genes in the set present in the analysis\n")
      cat("- **enrichmentScore**: raw enrichment score\n")
      cat("- **NES**: Normalized Enrichment Score — comparable across gene sets of different sizes\n")
      cat("- **pvalue**: nominal p-value\n")
      cat("- **p.adjust**: Benjamini-Hochberg adjusted p-value\n")
      cat("- **qvalue**: FDR q-value\n\n")
    }}
  }}
}}
```

```{{r gsea-table-{i}}}
# Display top GSEA results
if (!isTRUE(analysis_params$skip_gsea)) {{
  if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$gsea)) {{
    gsea_results <- as.data.frame(results[[{i}]]$gsea) %>%
        mutate(enrichmentScore=formatC(enrichmentScore, format="e", digits=2),
               NES=formatC(NES, format="e", digits=2),
               p.adjust=formatC(p.adjust, format="e", digits=2),
               qvalue=formatC(qvalue, format="e", digits=2))
    if (nrow(gsea_results) > 0) {{
      DT::datatable(gsea_results %>%
                   dplyr::select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue) %>%
                   head(20),
                   caption = "Top enriched gene sets")
    }} else {{
      cat("No significant enriched gene sets found")
    }}
  }} else {{
    cat("No GSEA results available for this comparison")
  }}
}}
```')
    
    rmd_content <- paste0(rmd_content, comparison_section)
  }
  
ipa_content <- '\n
## Utilizing IPA

The CSV Differential Abundance output from limma (available in the results directory
provided alongside this report *de_data/[comparison_name]_limma.csv*), can be
uploaded directly into **QIAGEN Ingenuity Pathway Analysis (IPA)** for self-exploration of
pathways predicted to be enriched by this experimental condition. Massey’s BISR provides
access to VCU’s license of IPA. If you do not already have an account associated with this
license, you may reach out to **morecockcm@vcu.edu** with your name, VCU health or VCU
email, and request for IPA. To perform a core abundance analysis, login with your
credentials here: **https://analysis.ingenuity.com/pa** and follow the instructions [here](https://qiagen.my.salesforce-sites.com/KnowledgeBase/KnowledgeNavigatorPage?id=kA41i000000L6rMCAS).

We host an annual hands-on training for IPA at the beginning of the fall semester. Please
email BISR if you would like to be a part of this training. In the meantime, QIAGEN has a
playlist of user-friendly tutorials available on Youtube titled “QIAGEN IPA Training
Videos” the **Qiagen Digital Insights Youtube** page.'

manuscript_content <- '\n
## Manuscript-Ready Text

### Methods
Prior to differential abundance analysis, missing intensity values were imputed
using the DEP R package [1] or an in-house method developed by the VCU Proteomics Core. Missing values in proteomics data-dependent
acquisition (DDA) experiments are commonly attributed to a missing-not-at-random
(MNAR) mechanism, where low-abundance peptides fall below the instrument
detection threshold. To address this, remaining missing values were imputed from
a down-shifted Gaussian distribution centered at the low end of the observed
intensity distribution.

Differential abundance analysis was carried out using limma v 1.44.0 [2].
Proteins were considered significantly differentially abundant (DAP) when they
exhibited an absolute log₂ fold-change ≥ 0.58 (1.5-fold) and an adjusted P-value
(FDR) < 0.05 (correction done with a Benjamini Hochberg). Volcano plots and
heatmaps were generated using the raw or log transformed intensity values and
visualized using R packages (ggplot2 and ggrepel), where the top up-regulated
and down-regulated proteins (based on adjusted P-value) were highlighted.
Heatmaps were produced using ComplexHeatmap after z-score transformation of the
filtered expression matrix. Gene Set Enrichment Analysis (GSEA) [3] for Gene
Ontology terms (GO) was performed using the clusterProfiler package [4] across
all proteins, regardless of significance. All computational analyses were
performed on VCU’s High Performance Research Computing cluster.

### References
1) Zhang X, Smits AH, van Tilburg GB, Ovaa H, Huber W, Vermeulen M. (2018).
Proteome-wide identification of ubiquitin interactions using UbIA-MS. Nature
Protocols 13, 530-550. https://doi.org/10.1038/nprot.2017.147

2) Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., Smyth, G.K., 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47. https://doi.org/10.1093/nar/gkv007

3) A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert, M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, & J.P. Mesirov, Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550, https://doi.org/10.1073/pnas.0506580102 (2005).

4) Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.

### Required Acknowledgements

Please include the following statements in your acknowledgements manuscript section:

- “Services in support of the research project were provided by the VCU Massey Comprehensive Cancer Center Bioinformatics Shared Resource. Massey is supported, in part, with funding from NIH-NCI Cancer Center Support Grant P30 CA016059.”

- “High Performance Computing resources provided by the High Performance Research Computing (HPRC) core facility at Virginia Commonwealth University (https://hprc.vcu.edu) were used for conducting the research reported in this work.”  
'
  # pasting together these last sections
  rmd_content <- paste0(rmd_content, ipa_content, manuscript_content)
  
  # Write the R Markdown file
  #fn <- glue('{report_prefix}_{timestamp}.Rmd')
  fn <- glue('{report_prefix}.Rmd')
  
  rmd_file <- file.path(output_dir, fn)
  writeLines(rmd_content[[1]], rmd_file)
  
  # Render the R Markdown to HTML
  rmarkdown::render(rmd_file, output_file = output_file, quiet = FALSE)
  
  # Return the path to the generated report
  return(output_file)
}