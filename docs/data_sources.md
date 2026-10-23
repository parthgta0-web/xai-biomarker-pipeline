# Data sources

The pipeline downloads public GEO series through `GEOparse` and caches them in `geo_cache/`. The cache is generated data and is intentionally ignored by Git.

| Preset | GEO series | Use in this project |
| --- | --- | --- |
| OSCC | GSE31056 | recurrence outcome |
| T2D | GSE22309, GSE25462 | skeletal-muscle diabetes comparison |
| Breast cancer | GSE27562, GSE47862 | PBMC classification |

Run `--inspect` before a dataset-specific analysis to review metadata and confirm that label rules still match GEO records.

## MitoCarta

The T2D mitochondrial subset expects a local MitoCarta 3.0 spreadsheet. The script reads worksheet **`A Human MitoCarta3.0`** and its **`Symbol`** column:

```bash
python xai_biomarker_pipeline.py --dataset t2d --mitocarta Human.MitoCarta3.0.xls
```

If a spreadsheet download is unavailable, a plain-text file with one gene symbol per line can be used instead.

## Entrez-to-symbol fallback

For platforms that do not supply usable gene symbols, the script downloads NCBI's public `Homo_sapiens.gene_info.gz` into `geo_cache/` and maps Entrez `GeneID` values to official symbols. A download failure stops with a clear message; manually placing the file in the cache permits the run to continue.

