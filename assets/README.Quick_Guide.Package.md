# Proteomics Pipeline: Quick Guide (R Package Install)

<img src="../assets/DAPRmd.png" alt="DAPRmd Logo" width="800"/><br>

---

## Overview

DAPRmd is now structured as an R package. Installing it is a two-step process:

1. **conda** — installs R and all R package dependencies (CRAN + Bioconductor)
2. **`devtools::install('.')`** — installs the DAPRmd package code itself

Conda is still required for R packages on RHEL 8 / glibc 2.28 systems (e.g. VCU Cardinal).
Precompiled CRAN/Bioconductor binaries from R's default repositories require glibc ≥ 2.32,
which RHEL 8 does not have. conda-forge builds its packages against glibc 2.17, so it works.

---

## Environment Setup

### What conda is still responsible for

- Providing a specific version of R (e.g. `r-base=4.5.3`)
- System libraries that compiled R packages depend on (e.g. `libcurl`, `libxml2`, `openssl`)
- **All R package dependencies** — required on RHEL 8 because precompiled CRAN/Bioconductor
  binaries need glibc ≥ 2.32, which the system does not provide

### What conda is no longer responsible for

- Installing the DAPRmd package itself — that is handled by `devtools::install('.')` in Step 2.

---

### Step 1: Create the conda environment with R and all dependencies

Use the provided environment file to create the full DAPRmd environment in one step:

```bash
conda env create -f envs/environment.cross_platform.yml
conda activate DAPRmd
```

> **Note:** The environment file installs R and all required CRAN and Bioconductor packages.
> For VCU Cardinal HPC, use `environment.vcu_cardinal.yml` instead for fully pinned builds.

---

### Step 2: Install the DAPRmd package itself

Once the conda environment is active, install just the DAPRmd package code:

```bash
# From the root of the proteomics_core_pipelines repo:
Rscript -e "devtools::install('.', dependencies=FALSE)"
```

The `dependencies=FALSE` flag is important — all R package dependencies were already
installed by conda in Step 1. This step only installs the DAPRmd package code itself.

---

### Step 3: Verify the installation

```bash
Rscript -e "library(DAPRmd); packageVersion('DAPRmd')"
```

---

## Alternative: Using a system module instead of conda for R

> **RHEL 8 / VCU Cardinal warning:** `module load R` gives you the system R, but precompiled
> CRAN/Bioconductor binaries from R's default repos require glibc ≥ 2.32. RHEL 8 ships
> glibc 2.28, so package installation will fail with a `GLIBC_2.32 not found` error.
> **The conda approach above is required on this system.**

On systems with glibc ≥ 2.32 (e.g. RHEL 9+), you can skip conda and use a module:

```bash
module load R/4.4.1
Rscript -e "install.packages('devtools')"
Rscript -e "devtools::install('.')"
```

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
Rscript de.regular.R \
  --counts human_counts.tsv \
  --samplesheet samplesheet.csv \
  --outdir human_results \
  --runid human_experiment \
  --annotation human
```

### Mouse Analysis

```bash
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
