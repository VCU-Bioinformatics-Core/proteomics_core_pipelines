# Proteomics Pipeline: Quick Guide

<img src="../assets/DAPRmd.png" alt="DAPRmd Logo" width="800"/><br>

---

## Overview

DAPRmd is structured as an R package. Installation involves two main steps:

1. **conda** creates an environment with a pinned R version and system libraries
2. **R package dependencies** are installed either automatically via `devtools` (Approach 1) or pre-installed from conda YAML files (Approach 2/3)

Three installation approaches are provided, in order of preference:

- **Approach 1** — conda creates a bare R environment; `devtools` resolves and installs all R dependencies at install time
- **Approach 2** — conda YAML file pre-installs R and all dependencies; then `devtools` installs DAPRmd code only
- **Approach 3** — same as 2 but dependencies are installed via explicit conda commands rather than a YAML file

When installing DAPRmd on the VCU Cardinal server (or similar systems), we recommend conda due to system compatibility constraints. Cardinal runs RHEL 8 with glibc 2.28, which is compatible with conda-forge packages (built against glibc 2.17) but not with standard precompiled CRAN/Bioconductor binaries. These require glibc ≥ 2.32 and will fail on Cardinal. Note: JP is looking into upgrading glibc on Cardinal.

---

## Environment & Package Setup 

| Section | Description | Link |
|---|---|---|
| Approach 1 | First try — minimal setup, devtools resolves all R dependencies automatically | [Approach 1](#approach-1-r-via-conda-with-devtools-dependency-installation-recommended) |
| Approach 2 | If Approach 1 fails or you want a fully reproducible environment with pinned package versions | [Approach 2](#approach-2-r-and-dependencies-via-conda-yaml-files-if-previous-fails) |
| Approach 3 | If Approach 2 fails or no YAML file is available for your system | [Approach 3](#approach-3-r-and-dependencies-via-conda-commands-if-previous-fails) |
| Verify the installation | Confirm DAPRmd installed correctly | [Verify](#verify-the-installation) |
| Why conda? | Explanation of why conda is used for package management | [Why conda?](#why-does-conda-handle-package-dependencies-in-approaches-23) |

### Approach 1: R via Conda with devtools Dependency Installation (Recommended)

Step 1: Create a bare conda environment with R.
```bash
conda create -n DAPRmd r-base=4.5.3 
```

Step 2: Activate the environment, clone the repo, and install DAPRmd (devtools will resolve and install all R dependencies automatically).
```bash
# activate the env
conda activate DAPRmd 

# clone the repo
git clone git@github.com:VCU-Bioinformatics-Core/proteomics_core_pipelines.git

# cd into the repo
cd proteomics_core_pipelines

# from the root of the proteomics_core_pipelines repo:
Rscript -e "devtools::install('.', dependencies=TRUE)"
```

### Approach 2: R and Dependencies via Conda YAML files (If previous fails)

Step 1: Create the conda environment from a YAML file — choose the one that matches your system.
```bash
# use on macOS or Linux with glibc ≥ 2.32 (e.g. RHEL 9+)
conda env create -f envs/environment.cross_platform.yml  

# use on VCU Cardinal or other RHEL 8 / glibc 2.28 systems
conda env create -f envs/environment.vcu_cardinal.yml    
```

Step 2: Activate the environment, clone the repo, and install DAPRmd (all R dependencies were already installed by conda in Step 1).
```bash
# activate the env
conda activate DAPRmd 

# clone the repo
git clone git@github.com:VCU-Bioinformatics-Core/proteomics_core_pipelines.git

# cd into the repo
cd proteomics_core_pipelines

# from the root of the proteomics_core_pipelines repo:
Rscript -e "devtools::install('.', dependencies=FALSE)"
```

### Approach 3: R and Dependencies via Conda Commands (If previous fails)

Step 1: Create a bare conda environment with R.
```bash
conda create -n DAPRmd r-base=4.5.3 
```

Step 2: Activate the environment
```bash
conda activate DAPRmd 
```

Step 3. Install CRAN/Bioconductor packages via conda

Each `conda install` command block below must be run separately. Installing all packages in a
single command can cause solver conflicts; running them one at a time ensures reliable
dependency resolution.

**CRAN packages** (via `conda-forge`):

```bash
conda install -c conda-forge \
  r-base=4.5.3 \
  r-circlize \
  r-data.table \
  r-devtools \
  r-dplyr \
  r-dt \
  r-futile.logger \
  r-ggplot2 \
  r-ggrepel \
  r-glue \
  r-haven \
  r-here \
  r-htmlwidgets \
  r-igraph \
  r-janitor \
  r-knitr \
  r-optparse \
  r-plotly \
  r-purrr \
  r-ragg \
  r-rcolorbrewer \
  r-readr \
  r-reticulate \
  r-rmarkdown \
  r-rvest \
  r-scales \
  r-stringr \
  r-textshaping \
  r-tibble \
  r-tidyr \
  r-tidyverse \
  r-tinytex \
  r-xml2
```

**Bioconductor packages requiring both `conda-forge` and `bioconda` channels** (run separately from the command above):

```bash
conda install -c conda-forge -c bioconda \
  bioconductor-annotationdbi \
  bioconductor-biomart \
  bioconductor-complexheatmap \
  bioconductor-go.db \
  bioconductor-gosemsim \
  bioconductor-limma \
  bioconductor-mzid \
  bioconductor-mzr \
  bioconductor-org.hs.eg.db \
  bioconductor-qfeatures \
  bioconductor-s4vectors \
  bioconductor-summarizedexperiment \
  bioconductor-vsn \
  bioconductor-affy
```

**Bioconductor packages from `bioconda` only** (run separately):

```bash
conda install -c bioconda \
  bioconductor-dose \
  bioconductor-msnbase \
  bioconductor-psmatch
```

**Bioconductor enrichment packages** (run separately — these depend on the packages above):

```bash
conda install -c bioconda \
  bioconductor-clusterprofiler \
  bioconductor-dep \
  bioconductor-enrichplot
```

Step 4: Clone the repo and install DAPRmd (all R dependencies were already installed by conda in Step 3)

```bash
# clone the repo
git clone git@github.com:VCU-Bioinformatics-Core/proteomics_core_pipelines.git

# cd into the repo
cd proteomics_core_pipelines

# from the root of the proteomics_core_pipelines repo:
Rscript -e "devtools::install('.', dependencies=FALSE)"
```

### Verify the installation

```bash
Rscript -e "library(DAPRmd); packageVersion('DAPRmd')"
```

### Why is conda integrated into this package, especially approaches 2/3?

- Provides a specific version of R (e.g. `r-base=4.5.3`)
- System libraries are easily installed (e.g. `libcurl`, `libxml2`, `openssl`)
- **R package dependencies** on CRAN/Bioconductor require glibc ≥ 2.32 but the VCU Cardinal server (and similar systems) only have access to older versions of glibc

---

## Preparing Your Data

Two input files are required in specific formats: the **Intensity Matrix** and the **Samplesheet**.

### 1. Intensity Matrix

**Required format:** Tab-Separated Values (`.tsv`)

This file contains the protein IDs and raw merged protein levels. The columns must be organized as follows:

- **PG.ProteinGroups:** The first column must contain protein identifiers in the form of UniProt/SwissProt IDs.
- **PG.Genes:** The second column must contain protein identifiers in the form of gene symbols.
- **Subsequent columns:** All remaining columns should contain raw intensity data for each sample. Column headers must exactly match the `SampleID` values in the samplesheet. Remove any additional columns (e.g. annotations, statistics) before running the pipeline.

| PG.ProteinGroups | PG.Genes | sample1_r1 | sample1_r2 | sample1_r3 | sample2_r1 | sample2_r2 | sample2_r3 | sample3_r1 | sample3_r2 | sample3_r3 | sample4_r1 | sample4_r2 | sample4_r3 |
| ---------------- | -------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| A0A087WV53       | SPEGNB   |    1243.58 |    2310.42 |    3428.67 |    4592.11 |    5635.74 |    6772.45 |    2384.91 |    3499.26 |    4621.08 |    5790.34 |    3435.89 |    4522.73 |
| A0A0B4J2D5       | P0DPI2   |    1189.24 |    2285.33 |    3375.19 |    4478.62 |    5592.47 |    6721.53 |    2290.16 |    3420.58 |    4518.47 |    5659.82 |    3392.41 |    4461.05 |
| A0A0B4J2F0       | PIGBOS1  |    1322.88 |    2399.75 |    3541.62 |    4684.91 |    5743.56 |    6833.10 |    2415.77 |    3556.28 |    4652.94 |    5734.62 |    3522.84 |    4613.25 |
| A0A0C4DH30       | IGHV3-16 |    1278.52 |    2366.42 |    3478.13 |    4555.08 |    5692.75 |    6811.44 |    2362.27 |    3518.32 |    4588.69 |    5688.11 |    3451.95 |    4582.63 |
| A0A0U1RRE5       | NBDY     |    1210.67 |    2333.91 |    3411.27 |    4518.45 |    5623.83 |    6745.92 |    2321.34 |    3465.28 |    4549.71 |    5679.53 |    3408.72 |    4528.84 |

### 2. Samplesheet

**Required format:** Comma-Separated Values (`.csv`)

This file contains metadata for each sample, including group identifiers and binary indicators for specific comparisons. The columns must be organized as follows:

- **SampleID:** Sample identifiers that exactly match the sample column headers in the intensity matrix.
- **GroupID:** Group identifiers for each sample (e.g., `experiment1`, `control1`). Used to define the experimental design in limma.
- **[comparison]:** Subsequent columns define pairwise comparisons. Column names follow the format `experiment_vs_control`. Use `1` for the experimental group, `0` for control, and leave blank to exclude a sample from that comparison.

| SampleID    | GroupID     | experiment1_vs_control1 | experiment2_vs_control2 |
| ----------- | ----------- | ----------------------- | ----------------------- |
| sample1_r1  | experiment1 | 1                       |                         |
| sample1_r2  | experiment1 | 1                       |                         |
| sample1_r3  | experiment1 | 1                       |                         |
| sample2_r1  | control1    | 0                       |                         |
| sample2_r2  | control1    | 0                       |                         |
| sample2_r3  | control1    | 0                       |                         |
| sample3_r1  | experiment2 |                         | 1                       |
| sample3_r2  | experiment2 |                         | 1                       |
| sample3_r3  | experiment2 |                         | 1                       |
| sample4_r1  | control2    |                         | 0                       |
| sample4_r2  | control2    |                         | 0                       |
| sample4_r3  | control2    |                         | 0                       |

---

## Running the Pipeline

### Human Analysis

```bash
Rscript inst/scripts/de.regular.R \
  --counts human_counts.tsv \
  --samplesheet samplesheet.csv \
  --outdir human_results \
  --runid human_experiment \
  --annotation human
```

### Mouse Analysis

```bash
Rscript inst/scripts/de.regular.R \
  --counts mouse_counts.tsv \
  --samplesheet samplesheet.csv \
  --outdir mouse_results \
  --runid mouse_experiment \
  --annotation mouse
```

### Example Run (Ovarian Cancer Dataset from MSstats)

The package includes a small example dataset derived from the MSstats ovarian cancer SRM dataset (14 proteins, 10 control samples, 6 tumor samples).

```bash
Rscript inst/scripts/de.regular.R \
  --counts inst/scripts/extdata/ovarian_intensity_matrix.csv \
  --samplesheet inst/extdata/ovarian_samplesheet.csv \
  --outdir ovarian_results \
  --runid ovarian_example \
  --annotation human \
  --imputation none \
  --skip-gsea
```

> `--skip-gsea` is used here because the dataset only contains 14 proteins, which is too few for meaningful Gene Ontology enrichment analysis.

### Arguments

| Flag | Default | Description |
|---|---|---|
| `-c, --counts` | — | Path to the intensity matrix **(Mandatory)** |
| `-s, --samplesheet` | — | Path to the samplesheet **(Mandatory)** |
| `-r, --runid` | — | Unique identifier for the run **(Mandatory)** |
| `-o, --outdir` | `./output` | Output directory |
| `-a, --annotation` | `mouse` | Genome annotation: `mouse` or `human` |
| `-i, --imputation` | `DEP-MinProb` | Imputation method: `DEP-MinProb`, `DEP-knn`, `DEP-bpca`, `DEP-QRILC`, `DEP-man`, custom (see `R/custom_imputation.R`), or `none` |
| `-q, --imputation-q` | `0.01` | Quantile cutoff for `DEP-MinProb` and `DEP-QRILC` |
| `--seed` | `42` | Random seed for stochastic imputation methods |
| `--heatmap-norm` | `zscore` | Heatmap row normalization: `zscore` or `none` |
| `--heatmap-top-n` | `1000` | Number of top proteins by CV in the global heatmap |
| `--gsea-ont` | `BP` | Gene Ontology category for GSEA: `BP`, `MF`, `CC`, or `ALL` |
| `--skip-gsea` | — | Skip GSEA analysis (flag) |
| `--skip-anova` | — | Skip one-way ANOVA (flag) |
| `--group-color1` | `#D55E00` | Hex color for the up-regulated group |
| `--group-color2` | `#0072B2` | Hex color for the down-regulated group |

---

## Pipeline Output

The pipeline generates an output directory (specified by `--outdir`) with the following structure:

```
[outDir]/
├── proteomics_analysis.html        ← final HTML report
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
