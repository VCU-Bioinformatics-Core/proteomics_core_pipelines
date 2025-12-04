# ==========================
# limma differential expression pipeline
# ==========================

library(here)
library(dplyr)
library(data.table)
library(tidyverse)
library(janitor)
library(scales)
library(clusterProfiler)
library(enrichplot)
library(readr)
library(DT)
library(limma)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(RColorBrewer)
library(purrr)
library(plotly)
library(stats)
library(orca)
library(reticulate)
library(optparse)
library(htmlwidgets)
library(ComplexHeatmap)
library(circlize)

# Connect to Ensembl (use the dataset for your species)
library(biomaRt)

# source helper functions
source(here::here('workflow/helpers.R'))
source(here::here("workflow/report_generator.phospho.R"))
#debug(generate_volcano)
#debug(run_analysis)
#debug(generate_heatmap)

# ==========================
# Command-line options
# ==========================
option_list = list(
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Required. Path to merged counts.tsv file"),
  make_option(c("-s", "--samplesheet"), type = "character", default = NULL,
              help = "Required. Path to samplesheet.csv file"),
  make_option(c("-o", "--outdir"), type = "character", default = "./output",
              help = "Output directory [default= %default]"),
  make_option(c("-r", "--runid"), type = "character", default = NULL,
              help = "Required. Unique run ID"),
  make_option(c("-a", "--annotation"), type = "character", default = "mouse",
              help = "Genome for annotation: 'mouse' or 'human' [default= %default]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

runID <- opt$runid
countData <- opt$counts
samplesheet <- opt$samplesheet
outDir <- opt$outdir
genome <- opt$annotation

# ==========================
# Debug mode
# ==========================
# debug <- TRUE
# if (debug){
#   runID <- "test_run"
#   countData <- "./example_counts.tsv"
#   samplesheet <- "./example_samplesheet.csv"
#   outDir <- "./test_output"
#   genome <- "mouse"
# }

debug <- TRUE
# debug <- FALSE
if (debug){
  # runID <- 'Test'
  # countData <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/raw/20251001_U266_UT_VS_S+B_LFQ_Report.csv'
  # samplesheet <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/raw/20251001_U266_UT_VS_S+B_LFQ_Report.samplesheet.csv'
  # outDir <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/analyses/'
  # annotation <- 'human'
  # genome <- 'human'
  
  runID <- 'Test'
  countData <- '/global/projects/proteomics_core/researcher_workspace/steven_grant/Bioinformatics/results/phospho_proteome/LFQ_phospho_4condition_DEV_Report_A.csv'
  samplesheet <- '/global/projects/proteomics_core/researcher_workspace/steven_grant/Bioinformatics/results/phospho_proteome/samplesheet.csv'
  outDir <- '/global/projects/proteomics_core/researcher_workspace/steven_grant/Bioinformatics/results/phospho_proteome/'
  annotation <- 'human'
  genome <- 'human'
}

# ==========================
# Load annotation DB
# ==========================
if (genome == "human") {
  if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
  annotation_db <- org.Hs.eg.db
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
} else if (genome == "mouse") {
  if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
  annotation_db <- org.Mm.eg.db
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
} else {
  stop("Invalid genome specified. Use 'mouse' or 'human'")
}

# ==========================
# Read and prepare data
# ==========================
out_dirs <- setup_directories(outDir)
full_prot_levels <- data.frame(read_csv(countData, col_names = TRUE))
full_prot_levels <- full_prot_levels %>% filter(!is.na(EG.PrecursorId)) # remove prots with missing names
print("WARNING JR: NEED TO ADDRESS THE DROPPING OF DUPLICATES")
full_prot_levels <- full_prot_levels %>% distinct(EG.PrecursorId, .keep_all = TRUE) # drop duplicate (temporarily)
colnames(full_prot_levels) <- sub("^X\\.\\d+\\.\\.", "", colnames(full_prot_levels)) # Remove the "X.#.." prefix only from columns that start with "X."
colnames(full_prot_levels) <- sub("\\.raw\\.EG\\.TotalQuantity\\..Settings.", "", colnames(full_prot_levels)) # Remove the "X.#.." prefix only from columns that start with "X."

# Extract a dataframe with just protein levels
# counts <- full_prot_levels %>% dplyr::select(where(is.numeric)) 
non_prot_cols <- c("EG.PrecursorId", 'PG.FASTAName')
counts <- full_prot_levels %>% dplyr::select(-any_of(non_prot_cols))

# Collect the comparisons from the samplesheet
comparisons_raw <- read.csv(samplesheet)
comparisons_cols <- colnames(comparisons_raw)[3:ncol(comparisons_raw)]
comparisons <- list()

for (comparisons_name in comparisons_cols) {

  comparisons_vector <- comparisons_raw[[comparisons_name]]
  
  exp_group <- unique(comparisons_raw$GroupID[comparisons_vector == 1 & !is.na(comparisons_vector)])
  ctrl_group <- unique(comparisons_raw$GroupID[comparisons_vector == 0 & !is.na(comparisons_vector)])
  
  if (length(exp_group) > 0 && length(ctrl_group) > 0) {
    comparisons[[length(comparisons) + 1]] <- list(
      name = comparisons_name,
      exp = exp_group,
      ctrl = ctrl_group
    )
  }
}

# Optional log2-transform (if values are not log-scaled)
intensity_matrix <- counts %>% mutate(across(where(is.numeric), ~ as.integer(round(.))))
intensity_matrix <- log2(intensity_matrix + 1)

# ==========================
# Set up the limma design matrix 
# ==========================
sample_info <- data.frame(sample = comparisons_raw$SampleID,
                          condition = comparisons_raw$GroupID)
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(factor(sample_info$condition))
limma_params <- list(E = intensity_matrix, design = design)

# ==========================
# Analysis loop
# ==========================

results <- list()
for (i in seq_along(comparisons)) {
  
  # run the analysis on the current samples
  curr_result <- run_analysis_phospho(comparisons[[i]], limma_params, intensity_matrix, out_dirs)
  
  # save the current results if successful
  if (!is.null(curr_result)){
    results[[i]] <- curr_result
  }
}

# ==========================
# Principal Component Analysis
# ==========================
print("# Principal Component Analysis")
input_pca_matrix <- as.matrix(intensity_matrix)

print("WARNING JR: THIS CODE IS LIKE A SCRAPPY IMPUTATOIN WITH MEAN VALUES")
input_pca_matrix <- apply(input_pca_matrix, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

print("WARNING JR: DIAGNOSE WHY WE STILL FIND NANS DESPITE THE PREVIOUS STEP")
if (anyNA(input_pca_matrix)){
  input_pca_matrix <- as.data.frame(input_pca_matrix) %>%
    dplyr::select(where(~ all(!is.na(.x) & !is.infinite(.x))))
}


# calculate PCA results and scores 
pca_results <- prcomp(input_pca_matrix, rank. = 3)
pc_scores <- pca_results$x

# calculate the variance explained from PCA
pca_var_pct <- round((pca_results$sdev^2)/sum(pca_results$sdev^2)*100, 1)

# create a dataframe to plot PC1 and PC2
two_comps_pca_df <- data.frame(Sample=rownames(pc_scores),
                     X=pc_scores[,1],
                     Y=pc_scores[,2],
                     Group=sample_info$condition)

# plot PC1 versus PC2 values
pca_plot <- ggplot(data=two_comps_pca_df, aes(x=X, y=Y, label=Sample, color=Group)) +
  geom_point(size = 1) +
  # geom_text_repel(aes(label = as.character(Sample)), show.legend = FALSE, size = 2.5) +
  #geom_text(aes(label = as.character(Sample)), show.legend = FALSE, size = 2.5) +
  xlab(paste("PC1 - ", pca_var_pct[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca_var_pct[2], "%", sep=""))
  #theme_bw()
print(pca_plot)
ggsave(pca_plot, filename = str_c(out_dirs$pca, "/PCA_plot.png"))

# # Interactive PCA 2D/3D # removing for proteomic core analyses
# print("# Interactive PCA 2D/3D")
# fig <- plot_ly(two_comps_pca_df, x = ~X, y = ~Y, color = ~Group,
#                colors = c('steelblue','firebrick','olivedrab', 'plum'),
#                type = 'scatter', mode = 'markers', size = 5.7)
# fig <- fig %>% layout(plot_bgcolor='#e5ecf6')
# export_plotly_to_html(fig, paste0(out_dirs$pca, "/allsamples_PCA_plot.html"))
# 
# components <- data.frame(pca_results$x)
# components$PC2 <- -components$PC2
# components$PC3 <- -components$PC3
# components <- cbind(components, Group=pca_df$Group)
# tot_explained_variance_ratio <- 100 * sum(summary(pca_results)[["importance"]]['Proportion of Variance',])
# fig3D <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group,
#                  colors = c('steelblue','firebrick','olivedrab', 'plum')) %>%
#   add_markers(size = 12) %>%
#   layout(title = paste('Total Explained Variance = ', tot_explained_variance_ratio),
#          scene = list(bgcolor = "#e5ecf6"))
# export_plotly_to_html(fig3D, paste0(out_dirs$pca, "/allsamples_PCA_plot3D.html"))

# ==========================
# Save RDS
# ==========================
print('Save RDS')
#rds <- list(results, comparisons, out_dirs, pca_plot, fig, fig3D)
rds <- list(results, comparisons, out_dirs, pca_plot)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#rds_path <- glue("analysis_results_{timestamp}.rds")
rds_path <- file.path(outDir, "data/analysis_results.rds")
saveRDS(rds, rds_path)

# ==========================
# Generate the HTML report
# ==========================
generate_report(rds_path, output_dir = outDir)
