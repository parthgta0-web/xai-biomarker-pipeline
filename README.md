# XAI Biomarker Pipeline

An explainable machine-learning workflow for discovering gene-expression biomarkers from public GEO microarray studies. The pipeline uses XGBoost and SHAP, with leakage-aware evaluation and differential-expression analysis.

It supports the three analyses developed in the accompanying Colab work:

- Oral squamous cell carcinoma (OSCC) recurrence — `GSE31056`
- Type 2 diabetes skeletal-muscle signatures — `GSE22309` and `GSE25462`, with an optional MitoCarta filter
- Breast-cancer PBMC classification — `GSE27562` and `GSE47862`

## Highlights

- Loads GEO studies and maps probes to gene symbols
- Merges studies, applies batch correction, and quantile-normalizes expression data
- Tunes an XGBoost classifier and reports hold-out metrics
- Ranks biomarkers with SHAP and compares all-gene and top-gene models
- Runs repeated cross-validation with SHAP feature selection performed inside each training fold
- Produces PCA, SHAP, ROC/confusion-matrix, metrics, and enrichment outputs

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\\Scripts\\Activate.ps1
pip install -r requirements.txt

# Offline smoke test with synthetic data
python src/xai_biomarker_pipeline.py --demo

# Inspect metadata before running an analysis
python src/xai_biomarker_pipeline.py --dataset oscc --inspect

# Run an analysis
python src/xai_biomarker_pipeline.py --dataset oscc
```

For the T2D MitoCarta analysis, supply a local MitoCarta 3.0 spreadsheet:

```bash
python src/xai_biomarker_pipeline.py --dataset t2d --mitocarta Human.MitoCarta3.0.xls
```

## Outputs

Each run writes a results folder containing hold-out and cross-validation metrics, ranked SHAP biomarkers, differential-expression statistics, optional enrichment results, and figures.

## Reproducibility notes

The workflow downloads public GEO data when run. It intentionally asks users to inspect sample metadata and verify class counts before interpreting a dataset-specific result. This is a research workflow, not a clinical diagnostic tool.

## Notebook

`notebooks/xai_biomarker_analysis.ipynb` is a cleaned, portable run notebook based on the original Colab workflow. The original public Colab remains available from the project link supplied with this work.

