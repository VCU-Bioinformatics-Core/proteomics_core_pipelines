#!/usr/bin/env Rscript
# msstats_to_daprmd.R
# Processes RawData.ovarian.csv (MSstats long format) through MSstats dataProcess(),
# then converts the protein-level output to DAPRmd input format:
#   - intensity matrix (wide): PG.ProteinAccessions, PG.Genes, one column per sample
#   - samplesheet: SampleID, GroupID, tumor_vs_control comparison column

library(MSstats)

# ── 1. Load raw data ──────────────────────────────────────────────────────────
# Run this script from the package root directory
extdata_dir <- "inst/scripts/extdata"
raw_path <- file.path(extdata_dir, "RawData.ovarian.csv")
raw <- read.csv(raw_path)
cat("Loaded", nrow(raw), "rows\n")

# ── 2. Summarize to protein level ─────────────────────────────────────────────
processed <- dataProcess(
  raw,
  normalization       = "equalizeMedians",
  summaryMethod       = "TMP",
  censoredInt         = "0",
  MBimpute            = TRUE,
  maxQuantileforCensored = 0.999
)

protein_data <- processed$ProteinLevelData
cat("Protein-level rows:", nrow(protein_data), "\n")
cat("Columns:", paste(colnames(protein_data), collapse = ", "), "\n")

# ── 3. Build wide intensity matrix ────────────────────────────────────────────
# MSstats reports log2 intensities in LogIntensities; DAPRmd expects raw values
# so we back-transform: 2^LogIntensities
# Sample ID = "Condition_BioReplicate" to match DAPRmd convention

protein_data$SampleID <- paste0(
  protein_data$GROUP, "_", protein_data$SUBJECT
)

wide <- reshape(
  protein_data[, c("Protein", "SampleID", "LogIntensities")],
  idvar     = "Protein",
  timevar   = "SampleID",
  direction = "wide"
)
colnames(wide) <- sub("^LogIntensities\\.", "", colnames(wide))

# Back-transform log2 intensities to raw intensities for DAPRmd
sample_cols <- setdiff(colnames(wide), "Protein")
wide[, sample_cols] <- 2^wide[, sample_cols]

# Add required DAPRmd metadata columns
wide <- data.frame(
  PG.ProteinAccessions = wide$Protein,
  PG.Genes             = wide$Protein,   # gene symbols not in this dataset
  wide[, sample_cols, drop = FALSE],
  check.names = FALSE
)

cat("\nIntensity matrix dimensions:", nrow(wide), "proteins x",
    length(sample_cols), "samples\n")
cat("Sample columns:", paste(sort(sample_cols), collapse = ", "), "\n")

# ── 4. Build samplesheet ──────────────────────────────────────────────────────
meta <- unique(protein_data[, c("SampleID", "GROUP")])
meta <- meta[order(meta$GROUP, meta$SampleID), ]
colnames(meta) <- c("SampleID", "GroupID")

# Create comparison: tumor (1) vs control (0)
meta$tumor_vs_control <- ifelse(meta$GroupID == "tumor", 1,
                         ifelse(meta$GroupID == "control", 0, NA))

cat("\nSamplesheet:\n")
print(meta)

# ── 5. Write outputs ──────────────────────────────────────────────────────────
out_dir <- extdata_dir

intensity_path  <- file.path(out_dir, "ovarian_intensity_matrix.csv")
samplesheet_path <- file.path(out_dir, "ovarian_samplesheet.csv")

write.csv(wide,  intensity_path,   row.names = FALSE)
write.csv(meta, samplesheet_path,  row.names = FALSE)

cat("\nWrote intensity matrix to:", intensity_path, "\n")
cat("Wrote samplesheet to:",      samplesheet_path, "\n")
