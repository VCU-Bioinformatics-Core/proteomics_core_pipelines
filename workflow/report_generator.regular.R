#!/usr/bin/env Rscript

# Required libraries
library(rmarkdown)
library(knitr)
library(dplyr)
library(glue)
library(limma)

# Function to generate automated R Markdown report
# source("report_generator.R"); generate_report('analysis.rds', output_dir = getwd()) # test below
generate_report <- function(analysis_results_path, output_dir = "./", report_prefix = "proteomics_analysis", analyst="Joaquin Reyna") {
  
  
  # Validate inputs
  if (!file.exists(analysis_results_path)) {
    stop("Analysis results file does not exist:", analysis_results_path)
  }

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
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
library(knitr)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(htmlwidgets)

# Load results
rds_data <- readRDS("{analysis_results_path}")
results <- rds_data[[1]]
comparisons <- rds_data[[2]]
out_dirs <- rds_data[[3]]
pca_plot <- rds_data[[4]]
#pca_plotly <- rds_data[[5]]
#pca_3d <- rds_data[[6]]
#annotation <- rds_data[[7]]
```

## Overview

**Analyst: {analyst}**

This report contains the results of differential abundance analysis for `r length(comparisons)` comparisons.
It starts with the **Sample Exploration using PCA** section that highlights sample similarity and
helps to explore and understand the relationships between samples. Next, we present the
**Differential Abundance Results** section with subsections for each comparison (more details
below). For further exploration of the results, please refer to the output directories
containing the raw data files.


## Pipeline

For this analysis we used the following steps:

1. **Differential abundance**: Limma was used to identify differentially abundant proteins.
2. **Functional analysis**: Gene Set Enrichment Analysis (GSEA) was performed using clusterProfiler.
3. **Thresholds**: Proteins with adjusted p-value < 0.05 and |FC| >= 1.5 (equivalent to log2(FC) = 0.58) were considered differentially abundant.

As part of this pipeline we produce the following files for your downstream use:

```
.
├── data
│   ├── de_data
│   │   └── limma_[comparison].csv
│   └── gsea_data
│       └── GO_Analysis_[comparison].csv
└── figures
    ├── gsea
    │   └── [comparison]_GSEA.png
    ├── heatmap
    │   └── [comparison]_heatmap.png
    ├── pca
    │   └── PCA_plot.png
    └── volcano
        └── [comparison]_volcano.png
```

## Sample Exploration using PCA

PCA was performed to visualize the overall patterns of protein levels across samples and
to identify potential batch effects or outliers. **What we expect:** Samples with similar
profiles should cluster together, while dissimilar samples are expected to
separate into distinct clusters.

```{{r pca-plot, fig.width=6, fig.height=4 }}
pca_plot
```

<!--Interactive PCA plots can also be found in the following files:

- 2D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot.html")`
- 3D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot3D.html")`-->

```{{r results-summary}}
# Create a summary table for all comparisons
summary_table <- data.frame(
  Comparison = character(),
  Experimental = character(),
  Control = character(),
  Total_DAPs = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  GSEA_Performed = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(comparisons)) {{
  if (!is.null(results[[i]])) {{
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
        GSEA_Performed = gsea_status,
        stringsAsFactors = FALSE
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
```\n')
  
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
  a **red horizontal line** to indicate the significant p-value threshold of \\< 0.05. Therefore,
  every point (protein) above that red line can be considered statistically signficant (note:
  before FDR correction). In addition, the **black vertical lines** indicate 1.5 fold change.
  Proteins highlighted in red are up-regulated in group1 compared to the group2. Proteins
  highlighted in blue are down-regulated in group1 compared to group2. Proteins in gray, do
  not meet the thresholds for both logFC and p-value.
- Data point: a protein
- X-axis: log2 fold-change (group1/group2)
- Y-axis: -log10(p-value)

```{{r volcano-{i}, out.width="80%", out.height="80%" }}
# Display volcano plot from file
volcano_path <- file.path(out_dirs$volcano, paste0("{name}_volcano.png"))
if (file.exists(volcano_path)) {{
  knitr::include_graphics(volcano_path)
}} else {{
  cat("Volcano plot not available for this comparison")
}}
```


**Heatmap**

- Description: a heatmap of **z-score normalized** intensity data with application of 
  **hierarchical clustering** by both samples and protein levels
- X-axis: samples
- Y-axis: proteins
- Color-scale: z-score normalized intensity levels 


```{{r heatmap-{i} }}
# Display heatmap from file
heatmap_path <- file.path(out_dirs$heatmap, paste0("{name}_heatmap.png"))
if (file.exists(heatmap_path)) {{
  knitr::include_graphics(heatmap_path)
}} else {{
  cat("Heatmap not available for this comparison")
}}
```


**Table of Differentially Abundant Proteins**

Table of the top differential abundant proteins with both nominal (**pvalue**)
and adjusted pvalues (**padj**).

```{{r top-daps-{i} }}
# Display top DAPs table

dap_flag = 0
if (!is.null(results[[2]]) && !is.null(results[[2]]$limma)) {{
  top_up <- results[[2]]$limma %>%
    filter(!is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC >= 0.58) %>%
    mutate(logFC = round(logFC, 2),
           pvalue = formatC(P.Value, format = "e", digits = 2),
           padj = formatC(adj.P.Val, format = "e", digits = 2)) %>%
    arrange(padj)
  
  top_down <- results[[2]]$limma %>%
    filter(!is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC <= -0.58) %>%
    mutate(logFC = round(logFC, 2),
           pvalue = formatC(P.Value, format = "e", digits = 2),
           padj = formatC(adj.P.Val, format = "e", digits = 2)) %>%
    arrange(padj)
  
  combined_sigs <- bind_rows(top_up, top_down)
  
  # update the flags
  if (nrow(combined_sigs) > 0) {{
    dap_flag <- 1
  }}

}}
```

```{{r daps-{i} }}

if (dap_flag > 0){{
  # DT::datatable(top_down %>% select(ENSEMBL_ID, SYMBOL, log2FoldChange, pvalue, padj, GENENAME),
  #             caption = "Top Down Regulated Proteins")
  
  DT::datatable(combined_sigs %>% dplyr::select(ensembl_gene_id, uniprotswissprot, logFC, pvalue, padj, uniprotswissprot),
                caption = "Differentially Abundance Proteins")
  
}} else {{
  cat("No differential abundance results available for this comparison\n")
}}
  
```

**Gene Set Enrichment Analysis**

- Description: Gene Set Enrichment Analysis using (GO) terms with both a bubble plot and
 table
- X-axis: the protein ratio (# genes related to Gene Set / total number of significant genes) to which
the term is enriched
- Y-axis: a given Gene Set

```{{r gsea-{i}, out.width="100%", out.height="100%"}}
# Display GSEA results from file
gsea_path <- file.path(out_dirs$gsea, paste0("{name}_GSEA.png"))
if (file.exists(gsea_path)) {{
  knitr::include_graphics(gsea_path)
}} else {{
  cat("GSEA results not available for this comparison")
}}
```


```{{r gsea-table-{i} }}
# Display top GSEA results
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

```')
    
    rmd_content <- paste0(rmd_content, comparison_section)
  }
  
ipa_content <- '\n
## Utilizing IPA

The CSV Differential Abundance output from limma (available in the results directory
provided alongside this report *./output/de_data/DESeq2_[comparison_name].csv*), can be
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

imputation_content <- '\n
## Imputation Reporting

The histograms below show the distribution of all log2 intensity values across
all proteins **before** (blue) and **after** (red) imputation. Imputed values
typically appear as a secondary peak shifted towards the lower end of the
distribution, reflecting the left-censored nature of missing-not-at-random
(MNAR) data in proteomics.

```{r imputation-histogram, fig.width=8, fig.height=5}
intensity_raw      <- rds_data[[5]]
intensity_imputed  <- rds_data[[6]]

raw_vals <- data.frame(
  value = unlist(intensity_raw),
  type  = "Pre-imputation"
)
raw_vals <- raw_vals[!is.na(raw_vals$value), ]

imp_vals <- data.frame(
  value = unlist(intensity_imputed),
  type  = "Post-imputation"
)

all_vals       <- rbind(raw_vals, imp_vals)
all_vals$type  <- factor(all_vals$type, levels = c("Pre-imputation", "Post-imputation"))

ggplot(all_vals, aes(x = value, fill = type)) +
  geom_histogram(alpha = 0.6, bins = 80, position = "identity") +
  scale_fill_manual(values = c("Pre-imputation" = "steelblue",
                               "Post-imputation" = "firebrick")) +
  labs(x = "log2 Intensity", y = "Count", fill = NULL,
       title = "Intensity distribution before and after imputation") +
  theme_bw() +
  theme(legend.position = "top")
```
'

manuscript_content <- '\n
## Manuscript-Ready Text

### Methods
Differential abundance analysis was carried out using limma v 1.44.0 [1].
Proteins were considered significantly differentially abundant (DAP) when they
exhibited an absolute log₂ fold-change ≥ 0.58 (1.5-fold) and an adjusted P-value
(FDR) < 0.05 (correction done with a Benjamini Hochberg). Volcano plots and
heatmaps were generated using the raw or log transformed intensity values and
visualized using R packages (ggplot2 and ggrepel), where the top up-regulated
and down-regulated proteins (based on adjusted P-value) were highlighted.
Heatmaps were produced using ComplexHeatmap after z-score transformation of the
filtered expression matrix. Gene Set Enrichment Analysis (GSEA) [2] for Gene
Ontology terms (GO) was performed using the clusterProfiler package [3] across
all proteins, regardless of significance. All computational analyses were
performed on VCU’s High Performance Research Computing cluster.

### References
1) Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., Smyth, G.K., 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47. https://doi.org/10.1093/nar/gkv007

2) A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert, M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, & J.P. Mesirov, Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550, https://doi.org/10.1073/pnas.0506580102 (2005).

3) Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.

### Required Acknowledgements

Please include the following statements in your acknowledgements manuscript section:

- “Services in support of the research project were provided by the VCU Massey Comprehensive Cancer Center Bioinformatics Shared Resource. Massey is supported, in part, with funding from NIH-NCI Cancer Center Support Grant P30 CA016059.”

- “High Performance Computing resources provided by the High Performance Research Computing (HPRC) core facility at Virginia Commonwealth University (https://hprc.vcu.edu) were used for conducting the research reported in this work.”  
'
  # pasting together these last sections
  # skipping IPA for now but will add on later is there is interest
  #rmd_content <- paste0(rmd_content, ipa_content, manuscript_content)
  rmd_content <- paste0(rmd_content, imputation_content, manuscript_content)
  
  # Write the R Markdown file
  #fn <- glue('{report_prefix}_{timestamp}.Rmd')
  fn <- glue('{report_prefix}.Rmd')
  
  rmd_file <- file.path(output_dir, fn)
  writeLines(rmd_content[[1]], rmd_file)
  
  # Render the R Markdown to HTML
  # rmarkdown::render(rmd_file, output_file = output_file, quiet = FALSE)
  
  # Return the path to the generated report
  return(output_file)
}

#debug = TRUE
#if (debug == TRUE){
#  rds_fn = '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/analyses/data/analysis_results.rds'
#  output_dir = '/global/home/reynaj/Projects/proteomics_core/analyst_workspace/pipeline_test_data/analyses/data/'
#  generate_report(rds_fn, output_dir = output_dir)
#}
