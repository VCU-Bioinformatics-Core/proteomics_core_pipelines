# Proteomics Pipeline: Quick Guide

<img src="../assets/DAPRmd.png" alt="DAPRmd Logo" width="800"/><br>

### Preparing Your R Environment

`renv` is used to create an isolated environment. If this is your first time working with renv then please start an R session and run the following:

```R
install.packages("renv")
```

Then go ahead and close this repository and cd inside
```bash
git clone https://github.com/VCU-Bioinformatics-Core/proteomics_core_pipelines
cd proteomics_core_pipelines
```

From there, start another R session and run the following:
```R
library(renv)
renv::restore()
```


### Preparing Your Data
Two input files are required in specific formats: the **Intensity Matrix** and the **Samplesheet**.

#### 1. Intensity Matrix
**Required format:** Tab-Separated Values (`.tsv`)
This file contains the protein ids and raw merged protein levels. The columns must be organized as follows:
- **PG.ProteinGroups:** The first column must contain protein identifiers in the form of UniProt/SwissProt IDs.
- **PG.Genes:** The second column must contain protein identifiers in the form of gene symbols.
- **Subsequent Columns:** All subsequent columns should contain the raw count data for each sample. The `header` values for these sample columns must match the sample identifiers (`SampleID`) used in the samplesheet. The pipeline expects integer counts, as is typical for RNA-seq data. While the script has the functionality to automatically parse and handle count data from various quantification tools, the users are still responsible for removing any additional columns besides the protein IDs and sample counts. It currently works best with the merged counts output from pipelines like [`nf-core's rnaseq Nextflow pipeline`](https://nf-co.re/rnaseq).

| PG.ProteinGroups | PG.Genes | sample1_r1 | sample1_r2 | sample1_r3 | sample2_r1 | sample2_r2 | sample2_r3 | sample3_r1 | sample3_r2 | sample3_r3 | sample4_r1 | sample4_r2 | sample4_r3 |
| ---------------- | -------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| A0A087WV53       | SPEGNB   |    1243.58 |    2310.42 |    3428.67 |    4592.11 |    5635.74 |    6772.45 |    2384.91 |    3499.26 |    4621.08 |    5790.34 |    3435.89 |    4522.73 |
| A0A0B4J2D5       | P0DPI2   |    1189.24 |    2285.33 |    3375.19 |    4478.62 |    5592.47 |    6721.53 |    2290.16 |    3420.58 |    4518.47 |    5659.82 |    3392.41 |    4461.05 |
| A0A0B4J2F0       | PIGBOS1  |    1322.88 |    2399.75 |    3541.62 |    4684.91 |    5743.56 |    6833.10 |    2415.77 |    3556.28 |    4652.94 |    5734.62 |    3522.84 |    4613.25 |
| A0A0C4DH30       | IGHV3-16 |    1278.52 |    2366.42 |    3478.13 |    4555.08 |    5692.75 |    6811.44 |    2362.27 |    3518.32 |    4588.69 |    5688.11 |    3451.95 |    4582.63 |
| A0A0U1RRE5       | NBDY     |    1210.67 |    2333.91 |    3411.27 |    4518.45 |    5623.83 |    6745.92 |    2321.34 |    3465.28 |    4549.71 |    5679.53 |    3408.72 |    4528.84 |

  
#### 2. Samplesheet:
**Required format:** Comma-Separated Values (`.csv`)
This file contains metadata for each sample, including group identifiers and binary indicators for specific comparisons. The columns must be organized as follows:
- **SampleID:** Sample identifiers that exactly match the sample column headers in the count matrix.
- **GroupID:** Group identifiers for each sample (e.g., `experiment1`, `control1`). These group identifiers are crucial for defining the experimental design in limma.
- [comparison]: Subsequent columns define pairwise comparisons. The column name should follow the format `experiment_vs_control`. For each comparison column, use `1` to indicate samples belonging to the experimental group, `0` for the control group, and leave the cell `blank` for samples to be excluded from that specific comparison. This design matrix setup allows the user to specify which samples are used for each comparison, providing flexibility in complex experimental designs.

| SampleID | GroupID      | experiment1_vs_control1	| experiment2_vs_control2	|
| -------- | ------------ | ------------------------- | ------------------------- |
| sample1_r1  | experiment1    | 1                  |              							|
| sample1_r2  | experiment1    | 1                  |							              |
| sample1_r3  | experiment1    | 1                  |							              |
| sample2_r1  | control1  | 0                 |                        	|
| sample2_r2  | control1  | 0                 |                        	|
| sample2_r3  | control1  | 0                 |                        	|
| sample3_r1  | experiment2 |                		    | 1                      	|
| sample3_r2  | experiment2 |                		    | 1                      	|
| sample3_r3  | experiment2 |                		    | 1                      	|
| sample4_r1  | control2 |                		| 0                      	|
| sample4_r2  | control2 |                		| 0                      	|
| sample4_r3  | control2 |                		| 0                      	|

### Running the Pipeline

The pipeline is executed using an R script. Example commands for running it on VCU HPRC (high performance research computing servers) are as follows:

#### Human Analysis

```bash
module load R/4.4.1

Rscript de.regular.R \
--counts human_counts.tsv \
--samplesheet samplesheet.csv \
--outdir human_results \
--runid human_experiment \
--annotation human
```

#### Mouse Analysis

```bash
module load R/4.4.1

Rscript de.regular.R \
--counts mouse_counts.tsv \
--samplesheet samplesheet.csv \
--outdir mouse_results \
--runid mouse_experiment \
--annotation mouse
```

#### Arguments
- `-c, --counts`: Path to the merged counts file **(Mandatory)**
- `-s, --samplesheet`: Path to the sample sheet file **(Mandatory)**
- `-o, --outdir`: Output directory (default: `./output`)
- `-r, --runid`: Unique identifier for the analysis run **(Mandatory)**
- `-a, --annotation`: Genome to use for annotation: `mouse` or `human` (default: `mouse`)
- `-i, --imputation`: Imputation method (default: `DEP-MinProb`). DEP methods: `DEP-MinProb`, `DEP-knn`, `DEP-bpca`, `DEP-QRILC`, `DEP-man`. Custom methods defined in `workflow/custom_imputation.R` (e.g. `3by3`). Use `none` to skip.
- `-q, --imputation-q`: Quantile cutoff for the left-censored distribution used by `DEP-MinProb` and `DEP-QRILC` (default: `0.01`)
- `--seed`: Random seed for reproducibility of stochastic imputation methods (default: `42`)
- `--heatmap-norm`: Heatmap row normalization — `zscore` (z-score normalize rows) or `none` (raw intensity) (default: `zscore`)
- `--heatmap-top-n`: Number of top proteins by CV to show in the global heatmap (default: `1000`)
- `--gsea-ont`: Gene Ontology category for GSEA: `BP`, `MF`, `CC`, or `ALL` (default: `BP`)
- `--skip-gsea`: Skip GSEA analysis for faster runs (flag, no value needed)
- `--skip-anova`: Skip one-way ANOVA for faster runs (flag, no value needed)
- `--group-color1`: Hex color for the up-regulated group in plots (default: `#D55E00`)
- `--group-color2`: Hex color for the down-regulated group in plots (default: `#0072B2`)

#### Current execution workflow (temporary)

At present, the pipeline is executed manually through an interactive RStudio session. This workflow exists to ensure a controlled runtime environment (via renv) and will likely be automated in the future.

The current steps are:

1. Launch **RStudio** on `cardinal.som.edu`.

2. Open the **proteomics_core_pipelines** GitHub repository as the active RStudio project.
This step is required to activate the project-specific **renv** environment and ensure all package versions are correctly resolved.

3. From within RStudio, open a terminal and navigate to the desired **analysis directory**.

4. Execute the pipeline using a wrapper script (example below).

5. Open the generated **R Markdown file** located in `<outdir>/*.Rmd` for inspection, debugging, or manual rendering.

##### Wrapper script invocation

The pipeline is currently launched via a shell wrapper that calls the R script with all required command-line arguments as show:

```
project_dir="<path-to-project-dir>"
Rscript <wrapper-script-path> \
  --counts <counts-path> \
  --samplesheet <samplesheet-path> \
  --outdir <output-dir> \
  --runid <arbitrary-id> \
  --annotation <human or mouse>
```

This approach allows the analysis to be parameterized from the command line while still relying on RStudio to manage the execution environment. The wrapper script is responsible only for argument passing; all analysis logic remains within the R pipeline itself. 

### Pipeline Output

The pipeline generates an output directory (specified by `--outdir`) with the following structure:

```
[outDir]/
├── proteomics_analysis.html        ← final HTML report (name controlled by --report-prefix)
├── data/
│   ├── de_data/
│   │   └── [comparison]_limma.csv
│   ├── gsea_data/
│   │   └── [comparison]_go_analysis_.csv
│   └── anova/
│       └── global_anova.csv
└── figures/
    ├── imputation/
    │   ├── global_imputation_histogram.png
    │   ├── global_imputation_distribution.png
    │   ├── global_imputation_missing_per_sample.png
    │   └── global_imputation_total_counts.png
    ├── pca/
    │   └── global_pca_plot.png
    ├── heatmap/
    │   ├── global_heatmap.png
    │   └── [comparison]_heatmap.png
    ├── volcano/
    │   └── [comparison]_volcano.png
    ├── ma/
    │   └── [comparison]_ma.png
    └── gsea/
        └── [comparison]_gsea.png
```

**Key output files:**
- `[comparison]_limma.csv` — full limma results (logFC, pvalue, padj, imputation_category per protein)
- `[comparison]_go_analysis_.csv` — GSEA Gene Ontology results (NES, pvalue, p.adjust, leading edge)
- `global_anova.csv` — one-way ANOVA results across all groups (generated when >2 groups)
- `global_pca_plot.png` — PCA of all samples colored by group
- `global_heatmap.png` — top proteins by coefficient of variation across all samples
- `global_imputation_*.png` — four diagnostic plots for the imputation step