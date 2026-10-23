# Limitations and interpretation

This repository is a student replication attempt, not validation of biomarkers or a clinical model.

1. **Cohorts differ from the papers.** The T2D replication uses two GEO series rather than the paper's three.
2. **Small samples create unstable estimates.** Scores can vary with split, seed, preprocessing, and class balance.
3. **Labels come from GEO metadata.** Explicit keyword rules require `--inspect`, but metadata can be ambiguous.
4. **No external validation cohort is used.** Internal scores do not establish clinical generalization.
5. **ComBat results may be optimistic.** The ComBat path applies label-aware correction before the split. Batch/label confounding can inflate its near-perfect T2D results; the MitoCarta z-score run is the main comparison.
6. **The T2D ComBat output has limited signal.** Only five genes have non-zero SHAP importance in the committed ComBat outputs; apparent top-20 filler genes and their p-values should not be over-interpreted.
7. **Breast-cancer batch handling can fall back to z-scoring.** ComBat may fail with missing values, and residual batch effects can remain; `ACTB` deserves caution as a possible batch-related feature.
8. **Feature selection is not biological validation.** SHAP ranks model contributions in this setup; it does not establish a causal mechanism or validated biomarker panel.

Use the committed CSVs for reproducibility, report preprocessing choices with metrics, and avoid cross-paper comparisons without accounting for cohort and protocol differences.

