# To-dos

Work to do on DAPRmd. Planned features (not yet designed or built) are described in
[Planned_Updates.md](Planned_Updates.md).

---

## 1. Test every protein-level run scenario

Run each scenario in [protein_run_scenarios.xlsx](protein_run_scenarios.xlsx) end to end
and record the result in its "Tested end-to-end" column. The phases come from the
spreadsheet's Phases sheet; scenario numbers match its `#` column.

For each scenario:

- Run the proteome pipeline with a config (simulated data is fine; reuse real data where
  available).
- Check that it finishes without errors and that the report renders, including the
  Statistical Model section.
- Where it applies, compare the results with hand-written limma code (same design,
  `eBayes(robust = TRUE)`).
- Update the spreadsheet: "Tested end-to-end" (Yes / No, with a note) and anything you
  found under "Supported".

### Done (rerun after every change as regression tests)

- [x] 1. One variable, pairwise (ovarian example)
- [x] 12. Two variables + covariates + blocking (simulated, matched hand-written limma)
- [x] 18. Numerical covariate (simulated age effect, matched hand-written limma)

### Phase 1: most common setups

- [ ] 4. One variable + covariates
- [ ] 5. Paired (patient in every group)
- [ ] 6. Repeated samples (patient in one group)
- [ ] 8. Two variables (pairs, pooled and interaction contrasts)
- [ ] 15. Old two-column samplesheet (`SampleID, GroupID`, no `primary_variables`)

### Phase 2: remaining combinations of the same building blocks

- [ ] 2. One variable, each vs control (`comparisons: vs_reference`)
- [ ] 3. One variable, pooled groups (`(drugA + drugB)/2 - control`)
- [ ] 7. One variable + covariates + blocking
- [ ] 9. Two variables + covariates
- [ ] 10. Two variables, patient in every combination
- [ ] 11. Two variables, patient in every level of one and a single level of the other

### Phase 3: edge cases and alternative entry points

- [ ] 13. Three or more variables
- [ ] 14. Unbalanced or missing combinations (check the warning for single-sample groups)
- [ ] 16. Old 0/1 comparison columns (deprecated path)
- [ ] 17. Command-line flags only, no config

### Phase 4: needs development first (not testable yet)

These aren't supported, so check only that each fails with a clear error until the
feature exists.

- [ ] 19. Two blocking variables
- [ ] 20. Additive model (see the formula entry in [Planned_Updates.md](Planned_Updates.md))
- [ ] 21. Subset of samples

### Notes

- Scenarios 5, 6 and 11 differ only in how patients relate to the groups. That changes
  the statistics, so test them separately even though the configs look alike.
- Known bug to avoid while testing PTM runs: `run_ptm_ida_pipeline` crashes with
  `ptm_imp_matrix not found` unless a custom imputation method is used (also present on
  `main`).
