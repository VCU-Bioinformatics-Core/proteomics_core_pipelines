# run from the package directory
Rscript inst/scripts/de.regular.R \
    --counts inst/scripts/extdata/ovarian_intensity_matrix.csv \
    --samplesheet inst/extdata/ovarian_samplesheet.csv \
    --outdir test/ovarian_results \
    --runid ovarian_example \
    --annotation human \
    --imputation none \
    --skip-gsea