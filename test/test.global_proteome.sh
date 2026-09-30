# run from the package directory
Rscript --vanilla inst/scripts/de.regular.R \
    --counts inst/scripts/extdata/ovarian_intensity_matrix.csv \
    --samplesheet inst/extdata/ovarian_samplesheet.csv \
    --outdir test/ovarian_results \
    --runid ovarian_example \
    --annotation human \
    --imputation none \
    --skip-gsea

# single-gene strip plot from the results above
Rscript --vanilla inst/scripts/plot_single_gene.R \
    --gene CD59 \
    --imputed test/ovarian_results/data/protein_imputed_matrix.csv \
    --limma test/ovarian_results/data/de_data/tumor_vs_control_limma.csv \
    --samplesheet inst/extdata/ovarian_samplesheet.csv \
    --raw test/ovarian_results/data/protein_raw_matrix.csv \
    --out test/ovarian_results/figures/CD59_strip_plot.png
