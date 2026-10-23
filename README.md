# XAI biomarker pipeline: replication attempt of three DTU papers

[![Python syntax check](https://github.com/parthgta0-web/xai-biomarker-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/parthgta0-web/xai-biomarker-pipeline/actions/workflows/ci.yml)

An exploratory student replication attempt using XGBoost and SHAP on public GEO microarray data. The pipeline uses XGBoost and SHAP, with leakage-aware evaluation and differential-expression analysis.

> **Publication status:** the original pipeline implementation was supplied by Yuvraj and then patched and run for this student project. Please keep the repository private until Yuvraj has confirmed permission to publish that implementation.

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

## Repository layout

```text
├── xai_biomarker_pipeline.py  # Main runnable pipeline
├── notebooks/   # Portable notebook entry point
├── results/     # Committed result tables and figures
├── docs/        # Limitations and data-source notes
├── requirements.txt
└── README.md
```

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\\Scripts\\Activate.ps1
pip install -r requirements.txt

# Offline smoke test with synthetic data
python xai_biomarker_pipeline.py --demo

# Inspect metadata before running an analysis
python xai_biomarker_pipeline.py --dataset oscc --inspect

# Run an analysis
python xai_biomarker_pipeline.py --dataset oscc
```

On Windows PowerShell, activate the virtual environment with:

```powershell
.venv\Scripts\Activate.ps1
```

For the T2D MitoCarta analysis, supply a local MitoCarta 3.0 spreadsheet:

```bash
python xai_biomarker_pipeline.py --dataset t2d --mitocarta Human.MitoCarta3.0.xls
```

## Outputs

Each run writes a results folder containing hold-out and cross-validation metrics, ranked SHAP biomarkers, differential-expression statistics, optional enrichment results, and figures.

The committed `results/` folders contain the supplied CSV summaries and figures. Downloaded GEO caches, archives, and local MitoCarta spreadsheets are excluded from version control.

## Data and study presets

The workflow obtains public expression data from NCBI GEO at runtime. Built-in presets cover:

| Preset | GEO series | Task |
| --- | --- | --- |
| `oscc` | GSE31056 | Oral squamous cell carcinoma recurrence |
| `t2d` | GSE22309, GSE25462 | Type 2 diabetes skeletal-muscle signature |
| `bc` | GSE27562, GSE47862 | Breast-cancer PBMC classification |

Use `--inspect` before each dataset-specific run to review sample metadata and confirm the expected classes.

## Reproducibility notes

The workflow downloads public GEO data when run. It intentionally asks users to inspect sample metadata and verify class counts before interpreting a dataset-specific result. This is a research workflow, not a clinical diagnostic tool.

## Notebook

`notebooks/xai_biomarker_analysis.ipynb` is a cleaned, portable run notebook based on the original Colab workflow. The original run history is available in [Google Colab](https://colab.research.google.com/drive/12mdYYT1KopItTEc3HFvGoFzzxvqPaVEh#scrollTo=qFlv5UN7fPVQ).

## Contributing

Contributions are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request.

## Reproducibility controls

The default run uses the full XGBoost grid and three repeated stratified 5-fold CV repeats. Use `--fast` for fixed reproducible settings (`n_estimators=100`, `max_depth=3`, `learning_rate=0.1`, `subsample=0.8`, and `colsample_bytree=0.8`) and `--cv-repeats N` to change the repeat count.

```bash
python xai_biomarker_pipeline.py --dataset bc --fast --cv-repeats 1
```

When a platform has Entrez IDs but no gene symbols, the script downloads the public NCBI `Homo_sapiens.gene_info.gz` file into `geo_cache/` and uses it as a fallback.

## Results at a glance

The following values are copied from the committed output CSV files. They are exploratory internal estimates, not clinical performance claims.

| Analysis | Hold-out accuracy (AUC) | CV accuracy (AUC) | Note |
| --- | --- | --- | --- |
| OSCC, GSE31056 | 0.750 (0.844) | 0.696 (0.755) | Paper: accuracy 0.70, AUC 0.89, 5-fold accuracy 0.766 |
| T2D, all genes + ComBat | 0.941 (1.000) | 0.969 (0.996) | Two datasets; not directly comparable to the paper cohort |
| T2D, MitoCarta + ComBat | 1.000 (1.000) | 0.988 (1.000) | Likely optimistic because of label leakage |
| T2D, MitoCarta + z-score | 0.941 (1.000) | 0.835 +/- 0.128 (0.837) | Main T2D comparison; paper accuracy 0.8055 |
| Breast cancer, all genes | 0.856 (0.920) | 0.888 +/- 0.044 (0.951) | Paper accuracy 0.9667 |
| Breast cancer, SHAP top-20 | 0.844 (0.892) | 0.861 +/- 0.033 (0.907) | Paper top-10 accuracy 0.9444 |

`results/oscc/`, `results/t2d_all_genes_combat/`, `results/t2d_mito_combat/`, `results/t2d_mito_zscore/`, and `results/bc/` contain the complete tables and figures. The T2D z-score top-20 model had CV accuracy 0.765 and AUC 0.820.

## Interpretation and gene overlap

- **OSCC:** leading SHAP genes (including `GPX7`, `PNPO`, `MRFAP1`, and `MED11`) did not overlap directly with the paper's `GRIA2`, `SERF2`, and `CARD8` set.
- **T2D:** `DNAJC4` and `BAX` occur in both the z-score result and the paper, but lower in the ranking; all-gene results also contain lower-ranked `ATP5B`, `ATP5C1`, and `ECHS1`.
- **Breast cancer:** `TMIGD2` is ranked first here and appears lower in the paper's top-10 list. `POMGNT2` appears in the paper; `ACTB` may reflect residual batch signal.

This project retains 482 shared MitoCarta genes, whereas the paper reports 409 out of 5,762 genes; the difference is expected because this replication uses two GEO datasets rather than three.

## Limitations

Small cohorts, metadata-derived labels, differing study composition, and no external validation limit interpretation. The ComBat path applies label-aware correction before splitting and can therefore inflate results when batch and label are entangled. The z-score T2D run is the more defensible comparison. See [limitations](docs/limitations.md) and [data sources](docs/data_sources.md) before reusing a result.

## Credits and references

The original pipeline implementation was supplied by **Yuvraj**. The student patches and documents that supplied work; they do not claim to have written the original pipeline from scratch.

- Dogra & Hasija (2025), DOI: [10.1016/j.imu.2025.102006](https://doi.org/10.1016/j.imu.2025.102006)
- Sharma, Saini & Hasija (2024), DOI: [10.3389/fgene.2024.1356746](https://doi.org/10.3389/fgene.2024.1356746)
- Kumar & Das (2023), DOI: [10.1016/j.compbiomed.2023.106861](https://doi.org/10.1016/j.compbiomed.2023.106861)
- Rath et al. (2021), MitoCarta3.0, DOI: [10.1093/nar/gkaa1011](https://doi.org/10.1093/nar/gkaa1011)

## Disclaimer

For educational and research exploration only. This project is not medical advice, a diagnostic tool, or evidence for clinical decision-making.

