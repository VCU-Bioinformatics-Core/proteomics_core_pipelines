# Proteomics Pipeline: Quick Guide

<img src="../assets/DAPRmd.png" alt="DAPRmd Logo" width="800"/><br>

---

## DAPRmd Environment Setup

This section describes how to set up a reproducible conda + R environment for running the
DAPRmd pipeline. R is compiled from source and installed directly into the conda environment
to avoid version conflicts with system R.

---

### 1. Create the conda environment

Create an empty conda environment named `DAPRmd`. No packages are installed at this stage —
R will be compiled from source in the next step.

```bash
conda create -n DAPRmd # creates an empty environment
conda activate DAPRmd
```

> **Note:** If you need system-level build dependencies (e.g. `gcc`, `make`, `readline-devel`),
> install them via `conda install -n DAPRmd <pkg>` before building R.

---

### 2. Install packages from a pre-built environment file (recommended)

If you are on a supported system, you can recreate the full DAPRmd environment in one step
using one of the provided environment files. This avoids the manual multi-step conda installs
described in the sections below.

Two environment files are available:

| File | When to use |
|------|-------------|
| `environment.cross_platform.yml` | General use — installs from `conda-forge` only, works on most Linux/macOS systems |
| `environment.vcu_cardinal.yml` | VCU Cardinal HPC — fully pinned build hashes for reproducibility on that cluster |

**Create the environment from a file:**

```bash
# Cross-platform (general use)
conda env create -f environment.cross_platform.yml

# VCU Cardinal HPC
conda env create -f environment.vcu_cardinal.yml
```

Then activate the environment:

```bash
conda activate DAPRmd
```

> **Note:** Both files create an environment named `DAPRmd`. If that name conflicts with an
> existing environment, rename it with `--name <new-name>` appended to the `conda env create`
> command.

If you cannot use an environment file (e.g. on a system with restricted network access or
non-standard architecture), follow the manual installation steps in sections 3–6 below.

---

### 3. Install CRAN packages via conda

Each `conda install` command below must be run separately. Installing all packages in a
single command can cause solver conflicts; running them one at a time ensures reliable
dependency resolution.

**CRAN packages** (via `conda-forge`):

```bash
conda install -c conda-forge \
  r-base=4.5.3 \
  r-biocmanager \
  r-data.table \
  r-dplyr \
  r-dt \
  r-ggplot2 \
  r-ggrepel \
  r-glue \
  r-here \
  r-htmlwidgets \
  r-janitor \
  r-knitr \
  r-optparse \
  r-plotly \
  r-purrr \
  r-rcolorbrewer \
  r-readr \
  r-reticulate \
  r-rmarkdown \
  r-scales \
  r-stringr \
  r-tibble \
  r-tidyr \
  r-tidyverse \
  r-tinytex \
  r-xml2 \
  r-haven \
  r-ragg \
  r-rvest \
  r-textshaping \
  r-devtools
```

**Bioconductor packages requiring both `conda-forge` and `bioconda` channels** (run separately from the command above):

```bash
conda install -c conda-forge -c bioconda \
  r-igraph \
  bioconductor-annotationdbi \
  bioconductor-biomart \
  bioconductor-gosemsim \
  bioconductor-mzr \
  bioconductor-vsn \
  bioconductor-affy \
  bioconductor-mzid \
  bioconductor-go.db \
  bioconductor-qfeatures \
  bioconductor-summarizedexperiment \
  bioconductor-org.hs.eg.db
```

**Bioconductor packages from `bioconda` only** (run separately):

```bash
conda install -c bioconda \
  bioconductor-psmatch \
  bioconductor-dose \
  bioconductor-msnbase
```

**Bioconductor enrichment packages** (run separately — these depend on the packages above):

```bash
conda install -c bioconda \
  bioconductor-enrichplot \
  bioconductor-clusterprofiler \
  bioconductor-dep
```

---

### Package dependency overview

| Layer | Packages |
|---|---|
| Report rendering | `rmarkdown`, `knitr` |
| Interactive output | `plotly`, `DT`, `htmlwidgets` |
| Visualization | `ggplot2`, `ggrepel`, `RColorBrewer`, `ComplexHeatmap` |
| Differential analysis | `limma`, `DEP` |
| Enrichment / annotation | `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`, `org.Mm.eg.db`, `AnnotationDbi`, `biomaRt` |
| Bioconductor infrastructure | `S4Vectors`, `SummarizedExperiment`, `BiocManager` |
| Data manipulation | `dplyr`, `tidyr`, `tibble`, `data.table`, `purrr`, `readr`, `janitor`, `scales` |
| Utilities | `glue`, `here`, `optparse`, `stringr`, `reticulate` |
| Reproducibility | `renv` |
| PDF output | `tinytex` |

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

The pipeline is executed via `Rscript`. Example commands for running on VCU HPRC are as follows:

### Human Analysis

```bash
module load R/4.4.1

Rscript de.regular.R \
  --counts human_counts.tsv \
  --samplesheet samplesheet.csv \
  --outdir human_results \
  --runid human_experiment \
  --annotation human
```

### Mouse Analysis

```bash
module load R/4.4.1

Rscript de.regular.R \
  --counts mouse_counts.tsv \
  --samplesheet samplesheet.csv \
  --outdir mouse_results \
  --runid mouse_experiment \
  --annotation mouse
```

### Arguments

| Flag | Default | Description |
|---|---|---|
| `-c, --counts` | — | Path to the intensity matrix **(Mandatory)** |
| `-s, --samplesheet` | — | Path to the samplesheet **(Mandatory)** |
| `-r, --runid` | — | Unique identifier for the run **(Mandatory)** |
| `-o, --outdir` | `./output` | Output directory |
| `-a, --annotation` | `mouse` | Genome annotation: `mouse` or `human` |
| `-i, --imputation` | `DEP-MinProb` | Imputation method: `DEP-MinProb`, `DEP-knn`, `DEP-bpca`, `DEP-QRILC`, `DEP-man`, custom (see `workflow/custom_imputation.R`), or `none` |
| `-q, --imputation-q` | `0.01` | Quantile cutoff for `DEP-MinProb` and `DEP-QRILC` |
| `--seed` | `42` | Random seed for stochastic imputation methods |
| `--heatmap-norm` | `zscore` | Heatmap row normalization: `zscore` or `none` |
| `--heatmap-top-n` | `1000` | Number of top proteins by CV in the global heatmap |
| `--gsea-ont` | `BP` | Gene Ontology category for GSEA: `BP`, `MF`, `CC`, or `ALL` |
| `--skip-gsea` | — | Skip GSEA analysis (flag) |
| `--skip-anova` | — | Skip one-way ANOVA (flag) |
| `--group-color1` | `#D55E00` | Hex color for the up-regulated group |
| `--group-color2` | `#0072B2` | Hex color for the down-regulated group |

### Current execution workflow (temporary)

At present, the pipeline is executed manually through an interactive RStudio session to ensure
a controlled runtime environment via `renv`. This will likely be automated in the future.

1. Launch **RStudio** on `cardinal.som.edu`.
2. Open the **proteomics_core_pipelines** repository as the active RStudio project to activate the `renv` environment.
3. From the RStudio terminal, navigate to the desired analysis directory.
4. Execute the pipeline via a wrapper script (example below).
5. Open the generated `.Rmd` file in `<outdir>/` for inspection or manual rendering.

#### Wrapper script invocation

```bash
Rscript <wrapper-script-path> \
  --counts <counts-path> \
  --samplesheet <samplesheet-path> \
  --outdir <output-dir> \
  --runid <arbitrary-id> \
  --annotation <human or mouse>
```

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
