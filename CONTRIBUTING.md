# Contributing

Thank you for helping improve the project.

## Before opening a pull request

1. Describe the purpose and scope of the change clearly.
2. Keep dataset-specific label logic explicit and document any changed assumptions.
3. Run the offline smoke test before submitting:

   ```bash
   python src/xai_biomarker_pipeline.py --demo
   ```

4. Do not commit downloaded GEO data, generated result directories, MitoCarta spreadsheets, or other large binary artifacts.

## Research reproducibility

This repository is intended for research and educational use. Changes that affect preprocessing, labels, validation, or feature selection should include enough context for reviewers to assess the effect on reproducibility and potential data leakage.

