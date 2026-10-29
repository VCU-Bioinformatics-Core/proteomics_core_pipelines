# Planned Updates

Changes we have discussed for DAPRmd but have **not** built yet. Each entry says what the
change would do, what it would affect in the code, and the order we plan to do it in.
Nothing here is available in the pipeline today.

---

## Passing a model formula through the config

**Status:** planned, not implemented (discussed 2026-10-29).

### Today's model

The pipeline always builds the same type of design in `build_design()`
([R/perform_limma_analysis.R](../R/perform_limma_analysis.R)):

```
~0 + GroupID + <categorical_covariates> + <numerical_covariates>
```

- `GroupID` gets one coefficient per combination of the `primary_variables` (a
  "cell-means" design). For example, with tissue (L, V) × response (NR, R) the groups are
  L_NR, L_R, V_NR and V_R.
- Covariates are adjusted for, never tested. Numerical covariates are centered on their
  mean and share one slope across all groups.
- Comparisons are written using group names (`V_R - V_NR`, or `(L_R + V_R)/2 - (L_NR + V_NR)/2`
  for a main effect).

The cell-means design gives the same fit as `~ tissue * response`; it's only labelled
differently. What it can't do is fit a **different** model, for example:

- additive (main effects only): `~ tissue + response`, with no interaction
- testing a numerical variable (proteins associated with age), rather than adjusting for it
- curved covariate effects: `splines::ns(age, df = 3)`
- group-specific slopes: `group:age`
- interactions between covariates: `batch:site`

### The idea

Add an optional config key that holds the model formula directly:

```yaml
design_formula: "~ tissue * response + batch + age"
```

The fit itself is easy: `model.matrix(as.formula(design_formula), data = samplesheet)`
and then `lmFit` as today. The difficulty is everything that assumes one coefficient per
group.

### What a free-form formula would affect

| Part of the pipeline | Today (cell-means) | With an arbitrary formula |
|---|---|---|
| Comparisons (`parse_comparisons`, `contrast_comparison` in [R/config.R](../R/config.R); `build_contrast_matrix` in [R/perform_limma_analysis.R](../R/perform_limma_analysis.R)) | written with group names, checked against the groups | must use coefficient names (`responseR`, `tissueV:responseR`) and be checked against `colnames(design)` |
| `comparisons: all` / `vs_reference` | generated from the groups | no longer mean anything |
| Volcano, heatmap and imputation categories | the samples on the + and − sides of a contrast are the "exp" and "ctrl" samples | an interaction or slope coefficient has no obvious set of samples to plot |
| Statistical Model report section ([R/report_model_section.R](../R/report_model_section.R)) | written specifically for cell-means plus covariates | needs a generic "here is your design and coefficients" version |
| Categorical vs numerical | set by `categorical_covariates` and `numerical_covariates` | the formula alone can't say whether `batch` is a category; it would need `factor(batch)` inside the formula, or keep the two lists |
| Config checks (`resolve_run_config` in [R/config.R](../R/config.R)) | unknown keys are rejected | the new key needs a matching pipeline argument in every pipeline that accepts it |
| Full-rank check, centering, `formula_label` / `numerical_covariate_means` attributes | done in `build_design()` | must be redone for a user formula |

**Easy to misread.** With `~ tissue * response`, the coefficient `responseR` is **not**
"responders vs non-responders". It is R vs NR **in the reference tissue only** (L, the
first level alphabetically). People often misread treatment-coded coefficients like this,
which is one reason the cell-means design is the default. A free-form formula exposes
every user to that mistake.

**Safety.** `as.formula()` on a config string can call any R function. That's acceptable
inside the core facility, but the pipeline should still check that the formula refers only
to samplesheet columns and a short list of allowed functions (`factor`, `log`, `ns`, ...).

### Three ways to offer it (safest first)

**1. Extra terms only.** Keep the cell-means part fixed and let the config add terms to it:

```yaml
extra_formula_terms: "splines::ns(age, df = 3) + batch:site"
```

The extra terms are adjusted for, like covariates, so comparisons, plots and the report
keep working. This covers curved covariate effects and interactions between covariates.
It's a small change.

**2. Named options for specific models.** Each common need gets its own key that the
pipeline fully understands:

```yaml
model: additive          # main effects only
test_numerical: [age]    # test the age slope (proteins associated with age)
ebayes_trend: true       # eBayes(trend = TRUE)
```

More work per feature, but each option gets proper plots, report text and checks.

**3. Expert mode with a full formula.** `design_formula` plus comparisons written over
coefficient names:

```yaml
design_formula: "~ tissue * response + batch"
comparisons:
  - name: R_vs_NR_in_L
    coef: responseR
  - name: interaction
    coef: "tissueV:responseR"
```

The pipeline would fit and test anything. Per-comparison heatmaps and imputation
categories would be skipped (or need extra config naming which groups to show), and the
report would only list the design and coefficients.

### Plan

1. **Option 1 first**: small, and it can't break existing runs.
2. **Option 2** for the specific models people actually ask for, starting with testing a
   numerical variable, then `ebayes_trend`, then `model: additive` (most useful when a
   combination of primary variables has no samples).
3. **Option 3 only if a power user really needs it.** Once a run uses a free-form
   formula, the reports can no longer be trusted to describe the model correctly.

### Before releasing any of these

- Existing configs must give byte-identical results (rerun the ovarian example and the
  simulated studies, scenarios 1, 12 and 18 in
  [protein_run_scenarios.xlsx](protein_run_scenarios.xlsx)).
- Check each new model against hand-written limma code, as was done for covariates and
  blocking.
- Update `README.md`, `assets/README.Quick_Guide.md` and the scenarios spreadsheet
  (row 20, "Additive model", is currently listed as not supported).

---

## limma variance settings in the config: `eBayes.trend`, `eBayes.robust`, `arrayWeights`

**Status:** planned, not implemented (requested 2026-10-29).

### Goal

Let users switch three limma variance settings from the run config. The defaults must
reproduce what the pipeline does today, so existing configs give identical results:

```yaml
eBayes.trend: false     # limma default; pipeline today: false
eBayes.robust: true     # limma default is false; pipeline today: true
arrayWeights: false     # pipeline today: not used
```

| Key | limma default | Pipeline today | Proposed default |
|---|---|---|---|
| `eBayes.trend` | `FALSE` | `FALSE` (not set) | `false` |
| `eBayes.robust` | `FALSE` | `TRUE` (hard-coded) | `true` |
| `arrayWeights` | not called | not called | `false` |

Note that `eBayes.robust` should default to **`true`**, the pipeline's current behaviour,
not limma's default of `FALSE`. Otherwise every existing run would change.

### What each setting does

limma's moderated t-test doesn't trust each protein's own variance on its own (with 3-4
samples per group it's very noisy). `eBayes` shrinks every protein's variance towards a
**prior** variance estimated from all proteins. The three settings change how that
borrowing works.

**`eBayes.trend`**: the prior variance depends on abundance.

- `false`: one prior variance for all proteins.
- `true`: the prior variance is a smooth curve of average log-intensity (`AveExpr`), so
  low-abundance proteins shrink towards a higher variance than high-abundance proteins.
- Why it matters for proteomics: low-intensity proteins are usually noisier (closer to
  the detection limit, more imputed values). With a single prior, these proteins get too
  little variance and too many false positives at low intensity, while high-intensity
  proteins are over-penalised. `trend = TRUE` is a common recommendation for proteomics
  (DEqMS takes the same idea further by modelling variance against the number of
  peptides/PSMs per protein).
- How to check it's needed: `plotSA(fit2)` plots residual standard deviation against
  average intensity. A clear slope means `trend = TRUE` is appropriate.

**`eBayes.robust`**: protects the prior from outlier proteins.

- `true` (pipeline today): proteins with unusually high or low variance (e.g. driven by a
  few odd samples or heavy imputation) are down-weighted when estimating the prior, and
  are themselves shrunk less, so they can't distort everyone else's variance.
- `false` (limma default): the classic method; slightly more powerful when there are no
  outlier variances.
- Combining it with `trend` is fine: `eBayes(fit, trend = TRUE, robust = TRUE)`.

**`arrayWeights`**: down-weights worse samples.

- Estimates one quality weight per sample (`arrayWeights(E, design)`) and passes it to
  `lmFit(..., weights = w)`. Noisy samples (poor digestion, low injection amount, a bad
  run) count less, instead of having to be dropped entirely.
- Most useful with a few clearly noisier samples; with uniformly good samples the weights
  are all close to 1 and nothing changes.
- The weights themselves are worth reporting: a sample with a very low weight is a
  quality flag in its own right.

### What it would affect in the code

`eBayes(..., robust = TRUE)` is currently hard-coded in three places, and all three must
read the new settings so the comparisons, the global F-test and the report agree:

| Where | What changes |
|---|---|
| `perform_limma_analysis()` in [R/perform_limma_analysis.R](../R/perform_limma_analysis.R) (per-comparison test) | `eBayes(fit2, trend = ..., robust = ...)` |
| `run_group_ftest()` in the same file (global moderated F-test) | same |
| `fit_limma_model()` in the same file | compute `arrayWeights` and pass `weights` to `lmFit` (see blocking note below) |
| `limma_params` in [R/run_proteome_da_pipeline.R](../R/run_proteome_da_pipeline.R) and [R/run_ptm_pipeline.R](../R/run_ptm_pipeline.R) | carry the three settings so the functions above can read them |
| Pipeline arguments | new arguments on both pipelines; `resolve_run_config` in [R/config.R](../R/config.R) rejects config keys that aren't pipeline arguments, so the keys only work once the arguments exist |
| CLI scripts in `inst/scripts/` | optional flags (`--ebayes-trend`, `--ebayes-robust`, `--array-weights`) |
| `run_config.yaml` | saved automatically like other settings, so every run records which settings it used |
| Parameters table in [R/generate_report_proteome.R](../R/generate_report_proteome.R) and [R/generate_report_ptm_std_ida.R](../R/generate_report_ptm_std_ida.R) | three new rows |
| Statistical Model section in [R/report_model_section.R](../R/report_model_section.R) | the prior-variance text currently says "robust empirical Bayes, `eBayes(robust = TRUE)`" and the equivalent-code block hard-codes `robust = TRUE`; both must follow the settings. With `arrayWeights`, add the weights to the model description and a table or plot of the per-sample weights |

**Config key names.** The dot in `eBayes.trend` is fine in YAML and in an R argument name.
The other config keys use snake_case (`blocking_variable`, `numerical_covariates`), so
`ebayes_trend`, `ebayes_robust` and `array_weights` would be more consistent. Decide
before implementing; the dotted form maps more directly onto limma's argument names.

**Not affected:** the PTM IDA pipeline ([R/run_ptm_ida_pipeline.R](../R/run_ptm_ida_pipeline.R))
uses per-site t-tests, not limma, so it should reject these keys with a clear error, the
same way it rejects blocking and covariates today.

### Interactions to handle

- **`arrayWeights` with a blocking variable.** `duplicateCorrelation` should be given the
  weights too, and ideally the weights and correlation are estimated in turn
  (weights → correlation → weights with the block). The simplest safe order:
  `w <- arrayWeights(E, design)`,
  `corfit <- duplicateCorrelation(E, design, block = block, weights = w)`,
  `lmFit(E, design, block = block, correlation = rho, weights = w)`.
- **`arrayWeights` with missing values.** It works with `NA`s, but with heavy
  missingness and no imputation the weights can be unstable. Log the weights and warn if
  any is extreme (e.g. below 0.2 or above 5).
- **`eBayes.trend` with imputation.** Imputed low values (e.g. DEP-MinProb) cluster at low
  intensity with artificially low variance, which can bend the trend. Worth checking with
  `plotSA` on a real dataset before recommending `trend = TRUE` by default.
- **`eBayes.trend` with few proteins.** The trend needs enough features to fit a curve;
  fine for a proteome, worth a check for small PTM datasets.

### Suggested extras (optional)

- Save the `plotSA` mean-variance plot in the report's QC section. It's the plot that tells
  users whether `eBayes.trend` is needed, so it helps with deciding on the setting.
- Report the estimated prior degrees of freedom (`fit2$df.prior`) and, with `robust = TRUE`,
  how many proteins were flagged as variance outliers.
- Allow `arrayWeights` options later (e.g. `var.design` to give the same weight to all
  samples from one batch), but start with the plain per-sample version.

### Plan

1. Add `eBayes.trend` and `eBayes.robust` together. They're the smallest change (only
   the `eBayes` calls and the report text), and the defaults keep results identical.
2. Add `plotSA` to the report, then run a few real datasets with `trend = TRUE` to decide
   whether it should become the default in a later release.
3. Add `arrayWeights`, including the blocking order above and a per-sample weight table in
   the report.

### Before releasing

- With default settings, rerun the regression scenarios (1, 12 and 18 in
  [protein_run_scenarios.xlsx](protein_run_scenarios.xlsx)) and confirm byte-identical
  output.
- For each non-default setting, compare against hand-written limma
  (`eBayes(..., trend = TRUE)`, `arrayWeights` + `lmFit(weights = )`) on a simulated
  dataset. For `arrayWeights`, simulate one deliberately noisy sample and check it gets a
  low weight.
- Update `README.md`, `assets/README.Quick_Guide.md` (settings reference), the commented
  example in `inst/extdata/ovarian_config.yaml`, and the "Options for every scenario" sheet
  of the scenarios spreadsheet.
