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
- **Approach 4** — uses renv which ensure compatibility with Rstudio instances (especially when a Docker/Singularity container is present)

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
conda install conda-forge::r-devtools
conda install -c conda-forge r-remotes

# clone the repo
git clone git@github.com:VCU-Bioinformatics-Core/proteomics_core_pipelines.git

# cd into the repo
cd proteomics_core_pipelines

# from the root of the proteomics_core_pipelines repo:
Rscript --vanilla -e "remotes::install_local('.', dependencies=TRUE)"

# try with devtools if the above fails
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

- **Approach 4** — use renv which ensure compatibility with Rstudio
This approach has to be applied to each project directory. 

Download the renv.lock file (from this repo) to your project directory
```
cd <your-project-dir>
wget https://raw.githubusercontent.com/VCU-Bioinformatics-Core/proteomics_core_pipelines/refs/heads/dev/renv.lock
```

Setup the environment and install the package (within an R session)
```
library(renv)
renv::init()
install.packages("remotes")
remotes::install_github("VCU-Bioinformatics-Core/proteomics_core_pipelines", ref="dev")
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

## Running an Analysis in 4 Steps

| Step | What you do | Section |
|---|---|---|
| 1 | Export the **intensity matrix** | [Intensity Matrix](#1-intensity-matrix) |
| 2 | Write a **samplesheet** (one row per sample: its ID and its variables) | [Samplesheet](#2-samplesheet) |
| 3 | Write a **run config** (YAML) that says which variables to compare and which groups to compare | [Run Config](#3-run-config) |
| 4 | Run one command | [Running the Pipeline](#running-the-pipeline) |

Most analyses only need the basics. If your samples have **more than one factor** (e.g. tissue × response), **several samples per patient**, or **batches**, see [Advanced Designs](#advanced-designs).

---

## Preparing Your Data

### 1. Intensity Matrix

**Required format:** Comma-Separated Values (`.csv`)

This file contains the protein IDs and raw merged protein levels. The columns must be organized as follows:

- **PG.ProteinAccessions:** The first column must contain protein identifiers in the form of UniProt/SwissProt IDs.
- **PG.Genes:** The second column must contain protein identifiers in the form of gene symbols.
- **Subsequent columns:** All remaining columns should contain raw intensity data for each sample. Column headers must exactly match the sample IDs in the samplesheet. Remove any additional columns (e.g. annotations, statistics) before running the pipeline.

| PG.ProteinAccessions | PG.Genes | sample1_r1 | sample1_r2 | sample1_r3 | sample2_r1 | sample2_r2 | sample2_r3 | sample3_r1 | sample3_r2 | sample3_r3 |
| -------------------- | -------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| A0A087WV53           | SPEGNB   |    1243.58 |    2310.42 |    3428.67 |    4592.11 |    5635.74 |    6772.45 |    2384.91 |    3499.26 |    4621.08 |
| A0A0B4J2D5           | P0DPI2   |    1189.24 |    2285.33 |    3375.19 |    4478.62 |    5592.47 |    6721.53 |    2290.16 |    3420.58 |    4518.47 |
| A0A0B4J2F0           | PIGBOS1  |    1322.88 |    2399.75 |    3541.62 |    4684.91 |    5743.56 |    6833.10 |    2415.77 |    3556.28 |    4652.94 |

### 2. Samplesheet

**Required format:** Comma-Separated Values (`.csv`)

The samplesheet describes **each sample**: one column of sample IDs, then one column per variable you know about the samples (treatment, tissue, patient, batch, ...).

- **Sample ID column** (e.g. `SampleID`): identifiers that exactly match the sample column headers in the intensity matrix.
- **Variable columns** (any names): e.g. `treatment`, `batch`.

| SampleID   | treatment | batch |
| ---------- | --------- | ----- |
| sample1_r1 | drugA     | b1    |
| sample1_r2 | drugA     | b2    |
| sample1_r3 | drugA     | b1    |
| sample2_r1 | drugB     | b2    |
| sample2_r2 | drugB     | b1    |
| sample2_r3 | drugB     | b2    |
| sample3_r1 | control   | b1    |
| sample3_r2 | control   | b2    |
| sample3_r3 | control   | b1    |

The samplesheet itself does not say what each column is for. The [run config](#3-run-config) does:

| Config key | What it is | Example |
|---|---|---|
| `sample_id` | the sample ID column | `SampleID` (the default) |
| `primary_variables` | the variable(s) you are comparing; their values are the **groups** | `[treatment]` → groups `drugA`, `drugB`, `control` |
| `covariates` | nuisance variables to adjust for (optional) | `[batch]` |
| `blocking_variable` | one variable whose samples are correlated, e.g. a patient ID (optional) | `PatientID` |

Columns not named in the config are ignored.

**Group names** (the values of the primary variables) must use only letters, digits, `.` or `_`, and must not start with a digit: `KO_1` and `h24` work; `KO-1` and `24h` do not. The pipeline stops with a clear message if a name is not allowed.

> **Older samplesheets** with a `GroupID` column still work: if `primary_variables` is not set, the `GroupID` column is used. Old `0`/`1` comparison columns in such samplesheets also still work when `comparisons` is not set, with a deprecation warning.

### 3. Run Config

**Required format:** YAML (`.yaml`)

One small text file per analysis holds every setting, so the whole analysis is described in one place and can be re-run exactly. Copy this template and edit it:

```yaml
# ---- inputs and outputs ----
run_id: my_experiment
protein_matrix_file: data/protein_intensities.csv
samplesheet_file: data/samplesheet.csv
out_dir: results/my_experiment

# ---- samplesheet columns ----
sample_id: SampleID               # column with the sample IDs
primary_variables: [treatment]    # what you are comparing; its values are the groups
covariates: [batch]               # optional: adjust for batch
# blocking_variable: PatientID    # optional: samples from the same patient

# ---- analysis settings (delete any line to use the default) ----
genome: human                     # human or mouse
imputation_method: DEP-MinProb    # see the settings table below
skip_gsea: false

# ---- which groups to compare ----
comparisons: all                  # all | vs_reference | a list (see below)
reference: control                # the control group
```

Relative paths are resolved from the folder you run the command in. All available settings are listed in the [Settings Reference](#settings-reference).

#### Choosing the comparisons

With `primary_variables: [treatment]` there are three groups (`drugA`, `drugB`, `control`):

| You write | Comparisons run |
|---|---|
| `comparisons: all` + `reference: control` | `drugA_vs_control`, `drugB_vs_control`, `drugB_vs_drugA` |
| `comparisons: vs_reference` + `reference: control` | `drugA_vs_control`, `drugB_vs_control` |
| a list of pairs (below) | exactly the pairs you list |

```yaml
comparisons:
  - [drugA, control]    # [experimental, control]
  - [drugB, control]
  - [drugA, drugB]
```

- **Direction:** in `[drugA, control]`, **positive logFC means higher in `drugA`**. Each comparison is named `<first>_vs_<second>`, and that name is used in the output file names.
- **`all` without `reference`:** the group that appears first in the samplesheet is used as the control.
- **`all` grows quickly:** *k* groups give *k(k−1)/2* comparisons (6 groups → 15), each with its own figures and GSEA. List the pairs you need when there are many groups.
- If `comparisons` is left out entirely, `all` is used.

---

## Running the Pipeline

### Try the example first (Ovarian Cancer Dataset from MSstats)

The package includes a small example dataset derived from the MSstats ovarian cancer SRM dataset (14 proteins, 10 control samples, 6 tumor samples) and a ready-made config, [`inst/extdata/ovarian_config.yaml`](../inst/extdata/ovarian_config.yaml). From the root of the repo:

```bash
Rscript inst/scripts/da.proteome.R --config inst/extdata/ovarian_config.yaml
```

The results are written to `test/ovarian_results/`. Its samplesheet has the columns `SampleID` and `condition`; the config sets `primary_variables: [condition]` and compares `tumor` against `control` and uses `skip_gsea: true` because 14 proteins are too few for meaningful Gene Ontology enrichment analysis.

### Run your own analysis

Keep the samplesheet and config in your project folder, then run the script for your data type:

```bash
cd <your-project-dir>

# protein-level differential abundance
Rscript <path-to-repo>/inst/scripts/da.proteome.R --config analysis.yaml

# PTM peptide-level differential abundance
Rscript <path-to-repo>/inst/scripts/da.ptm.R --config analysis.yaml

# PTM peptides adjusted for protein abundance (IDA)
Rscript <path-to-repo>/inst/scripts/ida.ptm.R --config analysis.yaml
```

For the PTM scripts, use `ptm_matrix_file` in the config instead of `protein_matrix_file`. `ida.ptm.R` needs both: `ptm_matrix_file` and `protein_matrix_file` (the imputed protein matrix).

### Changing a setting for one run

Any setting can be overridden on the command line without editing the config, e.g. a quick run without GSEA in a separate folder:

```bash
Rscript inst/scripts/da.proteome.R --config analysis.yaml --skip-gsea --outdir results_quick
```

You can also skip the config and pass everything as flags. Without a config, `comparisons` defaults to `all`, and the samplesheet needs `SampleID` and `GroupID` columns (`sample_id`, `primary_variables`, `covariates` and `blocking_variable` can only be set in a config):

```bash
Rscript inst/scripts/da.proteome.R \
  --counts mouse_counts.csv \
  --samplesheet samplesheet.csv \
  --outdir mouse_results \
  --runid mouse_experiment \
  --annotation mouse
```

### Re-running an analysis

Every run saves the settings it actually used, with the comparisons spelled out, to `<out_dir>/run_config.yaml`. To reproduce the run, pass that file back:

```bash
Rscript inst/scripts/da.proteome.R --config results/my_experiment/run_config.yaml
```

---

## Advanced Designs

The protein (`da.proteome.R`) and PTM (`da.ptm.R`) pipelines can also handle more complex studies. All samples are always analyzed together in one limma model, and each comparison is a question asked of that model. Sharing one model gives more statistical power than analyzing subsets separately.

| Your study has... | Add to the config | Section |
|---|---|---|
| more than one primary variable (e.g. tissue × response) | `primary_variables` | [Several primary variables](#several-primary-variables) |
| comparisons that pool or combine groups | formula comparisons | [Formula comparisons](#formula-comparisons) |
| several samples from the same patient | `blocking_variable` | [Several samples per patient](#several-samples-per-patient) |
| batches or other nuisance factors | `covariates` | [Batches and other covariates](#batches-and-other-covariates) |

### Worked example

Ventricular (`V`) and lumbar (`L`) samples from responders (`R`) and non-responders (`NR`). Each patient gave both a V and an L sample, and samples were run in two batches.

**Samplesheet:**

| SampleID | tissue | response | PatientID | batch |
| -------- | ------ | -------- | --------- | ----- |
| P1_V     | V      | R        | P1        | b1    |
| P1_L     | L      | R        | P1        | b1    |
| P2_V     | V      | R        | P2        | b2    |
| P2_L     | L      | R        | P2        | b2    |
| P5_V     | V      | NR       | P5        | b1    |
| P5_L     | L      | NR       | P5        | b1    |
| ...      | ...    | ...      | ...       | ...   |

**Run config:**

```yaml
run_id: csf_response
protein_matrix_file: data/protein_intensities.csv
samplesheet_file: data/samplesheet.csv
out_dir: results/csf_response
genome: human

sample_id: SampleID
primary_variables: [tissue, response]   # groups become V_R, L_R, V_NR, L_NR
covariates: [batch]                     # adjust for batch
blocking_variable: PatientID            # each patient gave several samples

comparisons:
  - [V_R, V_NR]                     # ventricular: responders vs non-responders
  - [L_R, L_NR]                     # lumbar: responders vs non-responders
  - [V_R, L_R]                      # responders: ventricular vs lumbar
  - [V_NR, L_NR]                    # non-responders: ventricular vs lumbar
  - name: R_vs_NR                   # all responders vs all non-responders
    contrast: "(V_R + L_R)/2 - (V_NR + L_NR)/2"
  - name: response_x_tissue         # is the response effect different in V and L?
    contrast: "(V_R - V_NR) - (L_R - L_NR)"
```

### Several primary variables

With one primary variable, its values are the groups. With several, `primary_variables: [tissue, response]` makes one group per combination, joined with `_`: `V_R`, `L_R`, `V_NR`, `L_NR`. Use these combined names in `comparisons`. The log lists every group with its number of samples.

### Formula comparisons

When a simple `[A, B]` pair is not enough, write a **formula** with the group names and give it a `name` (used in file names: letters, digits, `.`, `_`, `-`):

| Question | `contrast` |
|---|---|
| A vs B (same as the pair `[A, B]`) | `"A - B"` |
| the average of A and B vs C | `"(A + B)/2 - C"` |
| all responders vs all non-responders, each tissue counted equally | `"(V_R + L_R)/2 - (V_NR + L_NR)/2"` |
| does the R-vs-NR difference change between V and L? (interaction) | `"(V_R - V_NR) - (L_R - L_NR)"` |

- The groups with a `+` sign are the "experimental" side and the groups with a `-` sign are the "control" side. Positive logFC means higher on the `+` side, and the volcano, MA and heatmap figures use those samples.
- Use `/2`, `/3`, ... so each side is an **average**; otherwise the result is not a difference between groups (the pipeline warns about this).
- Typos in group names stop the run with a message listing the valid groups.

### Several samples per patient

If one patient (or animal, or pool) contributed more than one sample, set `blocking_variable` to the samplesheet column that identifies them:

```yaml
blocking_variable: PatientID
```

Samples from the same patient tend to look alike. Without `blocking_variable`, the pipeline would treat them as independent, and p-values for comparisons across patients (e.g. `R_vs_NR`) would be too small. With `blocking_variable`, limma estimates how strongly samples from the same patient are correlated (`duplicateCorrelation`) and accounts for it. The estimate is printed in the log and in the report; a value near 0 means blocking makes little difference.

Only one `blocking_variable` can be used.

### Batches and other covariates

```yaml
covariates: [batch]
```

Covariates adjust every comparison for nuisance factors such as batch or run date; the comparisons are written the same way. Covariates are treated as **categories** (batch `1`, `2`, `3` are three batches, not a number). Continuous covariates such as age are not supported yet.

If a covariate is mixed up with the groups (e.g. every responder was run in batch 1 and every non-responder in batch 2), the effects cannot be separated, and the pipeline stops and names the covariate. Remove it from `covariates` in that case.

### What changes in the output

- The report's **Analysis Parameters** table shows the limma design, the primary variables, the covariates and the blocking correlation, and each comparison section shows its limma contrast.
- The report's **Statistical Model** section (after the sample QC) writes out the model as an equation, defines every term for your study, shows the design matrix, each comparison as a formula and a contrast matrix, and the equivalent limma code.
- With more than two groups, the **global F-test** (`global_anova.csv`) lists proteins that differ between any of the groups. It uses the same model, so it also accounts for `covariates` and `blocking_variable`.

### PTM-IDA limitations

`ida.ptm.R` uses a different statistical test (per-site t-tests), so it accepts `primary_variables` and `[A, B]` pairs only. `covariates`, `blocking_variable` and formula comparisons stop the run with an error.

---

## Settings Reference

Each setting can go in the run config or be passed as a command-line flag. Command-line flags win.

| Config key | Flag | Default | Description |
|---|---|---|---|
| `run_id` | `-r, --runid` | — | Unique name for the run **(required)** |
| `protein_matrix_file` | `-c, --counts` | — | Protein intensity matrix **(required for `da.proteome.R`)**. For `ida.ptm.R`: the imputed protein matrix, flag `--protein-counts` |
| `ptm_matrix_file` | `-c, --counts` | — | PTM peptide intensity matrix **(required for `da.ptm.R` and `ida.ptm.R`)** |
| `samplesheet_file` | `-s, --samplesheet` | — | Samplesheet **(required)** |
| `out_dir` | `-o, --outdir` | `./output` | Output directory |
| `genome` | `-a, --annotation` | `mouse` | Genome annotation: `mouse` or `human` |
| `imputation_method` | `-i, --imputation` | `DEP-MinProb` | `DEP-MinProb`, `DEP-knn`, `DEP-bpca`, `DEP-QRILC`, `DEP-man`, custom (see `R/custom_imputation.R`), or `none` |
| `imputation_q` | `-q, --imputation-q` | `0.01` | Quantile cutoff for `DEP-MinProb` and `DEP-QRILC` |
| `imputation_seed` | `--seed` | `42` | Random seed for stochastic imputation methods |
| `heatmap_norm` | `--heatmap-norm` | `zscore` | Heatmap row normalization: `zscore` or `none` |
| `heatmap_top_n` | `--heatmap-top-n` | `1000` | Number of top proteins by CV in the global heatmap |
| `gsea_ont` | `--gsea-ont` | `BP` | Gene Ontology category for GSEA: `BP`, `MF`, `CC`, or `ALL` |
| `skip_gsea` | `--skip-gsea` | `false` | Skip GSEA analysis |
| `skip_anova` | `--skip-anova` | `false` | Skip the global F-test |
| `group_color1` | `--group-color1` | `#D55E00` | Hex color for the up-regulated group |
| `group_color2` | `--group-color2` | `#0072B2` | Hex color for the down-regulated group |
| `ptm_marker` | `--ptm-marker` | `Phospho \(STY\)` | PTM scripts only: pattern selecting PTM peptides |
| `comparisons` | — | `all` | Which groups to compare ([details](#choosing-the-comparisons)) |
| `reference` | — | — | The control group for `all` / `vs_reference` |
| `sample_id` | — | `SampleID` | Samplesheet column with the sample IDs |
| `primary_variables` | — | `GroupID` | Samplesheet column(s) being compared; their values are the groups ([details](#several-primary-variables)) |
| `covariates` | — | — | Columns to adjust for, e.g. batch ([details](#batches-and-other-covariates)) |
| `blocking_variable` | — | — | Column identifying samples from the same patient ([details](#several-samples-per-patient)) |
| — | `--config` | — | Path to the run config |

---

## Common Errors

The pipeline checks the samplesheet and config before the analysis starts, so most mistakes fail within seconds.

| Message starts with | What to do |
|---|---|
| `Group names (from primary_variables) must be valid R names` | Rename the groups, e.g. `24h` → `h24`, `KO-1` → `KO_1` |
| `Unknown key(s) in ...` | A config key is misspelled; the message lists the valid keys |
| `Missing required setting(s)` | Add the listed settings to the config or as flags |
| `Comparison [...] uses unknown group(s)` | Group names in `comparisons` must match the samplesheet groups exactly (or the combined names, with several `primary_variables`) |
| `Invalid contrast` | The formula uses a name that is not a group; the message lists the valid groups |
| `Column(s) not in the samplesheet` | A `sample_id`, `primary_variables`, `covariates` or `blocking_variable` name is misspelled; the message lists the samplesheet columns |
| `Column(s) assigned to more than one of ...` | A column can have only one role, e.g. not both a primary variable and a covariate |
| `Samplesheet SampleIDs and matrix columns do not match` | Make the sample IDs match the intensity matrix column headers |
| `The design matrix is not full rank` | A covariate is mixed up with the groups; remove it from `covariates` |
| `The PTM-IDA pipeline does not support ...` | Remove `covariates`, `blocking_variable` and formula comparisons for `ida.ptm.R` |
| Group names like `yes`, `no`, `on`, `off` behave oddly | YAML reads them as true/false; put them in quotes: `["yes", control]` |

---

## Pipeline Output

The pipeline generates an output directory (specified by `out_dir` / `--outdir`) with the following structure:

```
[outDir]/
├── proteomics_analysis.html        ← final HTML report (QC, statistical model, results per comparison)
├── run_config.yaml                 ← settings used for this run (re-run with --config)
├── data/
│   ├── de_data/
│   │   ├── [comparison]_limma.csv
│   │   └── master_query_table.csv
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

- `[comparison]_limma.csv` — full limma results (logFC, pvalue, padj, imputation_category per protein). Positive logFC = higher in the first group of the comparison.
- `master_query_table.csv` — for each protein, the comparisons in which it is significant
- `[comparison]_go_analysis_.csv` — GSEA Gene Ontology results (NES, pvalue, p.adjust, leading edge)
- `global_anova.csv` — global F-test across all groups, the limma version of a one-way ANOVA (generated when there are more than 2 groups)
- `run_config.yaml` — the settings used for this run
- `global_pca_plot.png` — PCA of all samples colored by group
- `global_heatmap.png` — top proteins by coefficient of variation across all samples
- `global_imputation_*.png` — four diagnostic plots for the imputation step
