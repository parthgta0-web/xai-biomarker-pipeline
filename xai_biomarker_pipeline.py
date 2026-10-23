#!/usr/bin/env python3
"""
XAI biomarker pipeline: XGBoost + SHAP on gene-expression microarray data.

Replicates the shared workflow of three DTU papers:
  1. Dogra & Hasija 2025      (OSCC recurrence, GSE31056)           -> --dataset oscc
  2. Sharma, Saini, Hasija 2024 (T2D skeletal muscle, mito genes)    -> --dataset t2d
  3. Kumar & Das 2023          (breast cancer PBMC, 2 GEO sets)      -> --dataset bc

Steps: load GEO -> probe-to-gene -> merge datasets -> batch correction ->
quantile normalisation -> (optional MitoCarta filter) -> 80/20 split ->
XGBoost + GridSearchCV -> metrics -> SHAP top-k genes -> retrain on top-k ->
leakage-free repeated CV -> differential expression (Welch t + BH FDR) ->
pathway enrichment (Enrichr via gseapy).

Usage
-----
  pip install xgboost shap scikit-learn statsmodels scipy matplotlib pandas GEOparse gseapy combat openpyxl xlrd
  python xai_biomarker_pipeline.py --demo                     # sanity check, no internet needed
  python xai_biomarker_pipeline.py --dataset oscc --inspect   # look at GEO metadata FIRST
  python xai_biomarker_pipeline.py --dataset oscc
  python xai_biomarker_pipeline.py --dataset t2d --mitocarta Human.MitoCarta3.0.xls
  python xai_biomarker_pipeline.py --dataset bc
  python xai_biomarker_pipeline.py --gse GSE12345 --pos "tumor" --neg "normal"   # any other dataset

IMPORTANT: sample labels are parsed from GEO metadata with keyword matching.
Always run with --inspect first and check the class counts it prints against
the paper (oscc: 38 vs 58, t2d: 70 vs 107, bc: 252 vs 194). Adjust --pos/--neg/--drop if needed.
"""
import argparse
import gzip
import re
import urllib.error
import urllib.request
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.metrics import (accuracy_score, confusion_matrix, f1_score,
                             precision_score, recall_score, roc_auc_score, roc_curve)
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, train_test_split
from statsmodels.stats.multitest import multipletests
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")

# label = 1 for "pos" (disease / recurrent / cancer), 0 for "neg". Negatives are checked first
# so that e.g. "non-recurrent" is not caught by "recurr". "drop" removes samples entirely.
PRESETS = {
    # outcome codes: DOD/AWD = disease came back (38); ANED/DOC = no recurrence (58)
    "oscc": dict(gse=["GSE31056"], label_col="outcome",
                 pos=["dod", "awd"], neg=["aned", "doc"],
                 drop=[], expected="38 recurrent vs 58 non-recurrent/normal"),
    # Per-dataset rules read from the real GEO metadata:
    #  GSE22309: use only samples BEFORE insulin; diabetic = 1, insulin sensitive = 0 (insulin resistant dropped)
    #  GSE25462: family-history column; "DM" = 1 (10), all other rows = 0 (40)
    #  GSE18732: glycemiagroup 3 = T2D, 1 = normal, 2 = IGT dropped (platform needs gene mapping, so off by default)
    "t2d": dict(gse=["GSE22309", "GSE25462"],
                per_gse={
                    "GSE22309": dict(label_col="source_name", pos=["diabetic individuals before"],
                                     neg=["insulin sensitive individuals before"], drop=["after"]),
                    "GSE25462": dict(label_col="family history", pos=["dm"], neg=["family history"], drop=[]),
                    "GSE18732": dict(label_col="glycemiagroup", pos=["3"], neg=["1"], drop=["2"]),
                },
                pos=[], neg=[], drop=[], expected="2 datasets: 25 T2D vs 60 control (paper total with GSE18732: 70 vs 107)"),
    "bc": dict(gse=["GSE27562", "GSE47862"],
               per_gse={
                   "GSE27562": dict(label_col="phenotype", pos=["malignant", "benign", "pre-surgery"],
                                      neg=["normal"], drop=["ectopic", "post-surgery"]),
                   "GSE47862": dict(label_col="developed_breast_cancer", pos=["yes"], neg=["no"], drop=[]),
               },
               pos=[], neg=[], drop=[], expected="252 breast cancer(+benign) vs 194 normal"),
}

GRID = dict(n_estimators=[50, 100, 200], max_depth=[3, 5], learning_rate=[0.05, 0.1, 0.2],
            subsample=[0.8], colsample_bytree=[0.8])
FAST_PARAMS = dict(n_estimators=100, max_depth=3, learning_rate=0.1,
                   subsample=0.8, colsample_bytree=0.8)
GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"


# ----------------------------------------------------------------------------- data loading
def _text_cols(pheno):
    keep = [c for c in pheno.columns
            if c in ("title", "source_name_ch1", "description") or c.startswith("characteristics")]
    return pheno[keep].astype(str)


def label_samples(pheno, pos, neg, drop, label_col=None):
    if label_col:   # read ONLY the metadata column(s) whose name contains label_col
        cols = [c for c in pheno.columns if label_col.lower() in c.lower()]
        if not cols:
            raise SystemExit(f"No metadata column contains '{label_col}'. Columns: {list(pheno.columns)}")
        text = pheno[cols].astype(str).apply(lambda r: " | ".join(r.values).lower(), axis=1)
    else:
        text = _text_cols(pheno).apply(lambda r: " | ".join(r.values).lower(), axis=1)

    def has(t, keys):
        return any(re.search(r"\b" + re.escape(k), t) for k in keys)

    def f(t):
        if has(t, drop):
            return np.nan
        if has(t, neg):
            return 0
        if has(t, pos):
            return 1
        return np.nan

    return text.map(f)


def clean_symbol(s):
    s = str(s)
    if " // " in s and "///" not in s:           # Affy ST style "NM_x // SYMBOL // desc"
        s = s.split(" // ")[1]
    s = re.split(r"\s*///\s*|\s*//\s*", s)[0].strip()
    s = re.sub(r"[\[\]<>]", "_", s)               # xgboost forbids these in feature names
    return "" if s.lower() in ("nan", "---", "none") else s


def load_ncbi_gene_info(cache):
    """Return an Entrez GeneID -> official symbol mapping, downloading the NCBI file if needed."""
    path = Path(cache) / "Homo_sapiens.gene_info.gz"
    if not path.exists():
        print("Downloading NCBI Homo_sapiens.gene_info.gz for Entrez-to-symbol mapping...")
        try:
            urllib.request.urlretrieve(GENE_INFO_URL, path)
        except (urllib.error.URLError, OSError) as exc:
            raise SystemExit(
                "Could not download NCBI gene_info. Check your connection or manually place "
                f"Homo_sapiens.gene_info.gz in {Path(cache).resolve()}: {exc}"
            ) from exc
    try:
        info = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, usecols=["GeneID", "Symbol"])
    except (OSError, ValueError, gzip.BadGzipFile) as exc:
        raise SystemExit(f"Could not read {path}: {exc}") from exc
    return info.dropna().drop_duplicates("GeneID").set_index("GeneID")["Symbol"]


def load_gse(gse_id, cache, inspect=False):
    import GEOparse
    gse = GEOparse.get_GEO(geo=gse_id, destdir=str(cache), silent=True)
    pheno = gse.phenotype_data
    if inspect:
        print(f"\n=== {gse_id}: {len(pheno)} samples, platform(s): {list(gse.gpls)} ===")
        for c in _text_cols(pheno).columns:
            print(f"-- {c}")
            print(pheno[c].value_counts().head(8).to_string())
        return None, pheno
    expr = gse.pivot_samples("VALUE")             # probes x samples
    ann = list(gse.gpls.values())[0].table
    cands = ["gene symbol", "gene_symbol", "symbol", "ilmn_gene", "gene_assignment", "gene"]
    col = next((c for c in ann.columns if c.lower() in cands), None)
    id_col = "ID" if "ID" in ann.columns else ann.columns[0]
    if col is None:
        entrez_col = next((c for c in ann.columns if c.lower().replace("_", " ") in
                           ("entrez gene id", "entrez geneid", "gene id")), None)
        if entrez_col is None:
            raise SystemExit(f"{gse_id}: no gene-symbol or Entrez-ID column in platform table. "
                             f"Columns: {list(ann.columns)}")
        print(f"{gse_id}: platform has no gene symbols; mapping Entrez IDs through NCBI gene_info.")
        entrez = ann.set_index(id_col)[entrez_col].astype(str).str.split(r"\\s*///\\s*|\\s*//\\s*").str[0]
        mapping = entrez.map(load_ncbi_gene_info(cache)).map(clean_symbol)
    else:
        mapping = ann.set_index(id_col)[col].map(clean_symbol)
    mapping.index = mapping.index.astype(str)
    expr.index = expr.index.astype(str)
    mapping = mapping.reindex(expr.index)
    ok = mapping.fillna("") != ""
    genes = expr[ok].groupby(mapping[ok]).mean()  # collapse probes -> genes (mean)
    if genes.max().max() > 50:                    # looks un-logged
        genes = np.log2(genes.clip(lower=1))
    return genes.T, pheno                         # samples x genes


def quantile_normalize(X):
    A = X.values
    order = np.argsort(A, axis=1)
    mean_sorted = np.take_along_axis(A, order, axis=1).mean(axis=0)
    out = np.empty_like(A)
    np.put_along_axis(out, order, np.broadcast_to(mean_sorted, A.shape), axis=1)
    return pd.DataFrame(out, index=X.index, columns=X.columns)


def batch_correct(X, batch, y, method):
    if batch.nunique() < 2 or method == "none":
        return X
    if method == "combat":
        try:
            from combat.pycombat import pycombat
            out = pycombat(X.T, batch.tolist(), mod=[y.tolist()])   # y kept as covariate
            return out.T
        except Exception as e:
            print(f"[warn] ComBat unavailable/failed ({e}); falling back to per-batch z-score")
    # location-scale: standardise each gene within each batch (what Kumar & Das describe)
    Z = X.copy()
    for b in batch.unique():
        idx = batch == b
        Z.loc[idx] = (X.loc[idx] - X.loc[idx].mean()) / X.loc[idx].std(ddof=0).replace(0, 1)
    return Z


def build_dataset(gse_ids, pos, neg, drop, cache, batch_method, mitocarta=None, label_col=None, per_gse=None):
    Xs, ys, bs = [], [], []
    for g in gse_ids:
        X, pheno = load_gse(g, cache)
        r = (per_gse or {}).get(g, dict(label_col=label_col, pos=pos, neg=neg, drop=drop))
        y = label_samples(pheno, r["pos"], r["neg"], r["drop"], r["label_col"]).reindex(X.index)
        print(f"{g}: {X.shape[1]} genes, labels -> {y.value_counts(dropna=False).to_dict()}")
        keep = y.notna()
        Xs.append(X[keep]); ys.append(y[keep].astype(int)); bs.append(pd.Series(g, index=X.index[keep]))
    X = pd.concat(Xs, join="inner")               # common genes only
    y, batch = pd.concat(ys), pd.concat(bs)
    print(f"Common genes across datasets: {X.shape[1]}")
    X = batch_correct(X, batch, y, batch_method)
    X = quantile_normalize(X)
    if mitocarta:
        mito = load_gene_list(mitocarta)
        X = X[[g for g in X.columns if g in mito]]
        print(f"After MitoCarta filter: {X.shape[1]} genes (paper 2: 409)")
    return X, y, batch


def load_gene_list(path):
    p = Path(path)
    if p.suffix.lower() in (".xls", ".xlsx"):
        df = pd.read_excel(p, sheet_name="A Human MitoCarta3.0")
        return set(df["Symbol"].astype(str))
    return set(l.strip() for l in p.read_text().splitlines() if l.strip())


def demo_dataset(seed=0, n=140, p=600, k_true=10):
    rng = np.random.default_rng(seed)
    y = pd.Series((rng.random(n) < 0.4).astype(int))
    batch = pd.Series(np.where(np.arange(n) < n // 2, "demoA", "demoB"))
    X = rng.normal(8, 1, size=(n, p))
    X[:, :k_true] += np.where(rng.random(k_true) < 0.5, 1, -1) * 1.3 * y.values[:, None]
    X[batch.values == "demoB"] += 0.8                     # fake batch effect
    cols = [f"GENE{i}" for i in range(p)]
    idx = [f"S{i}" for i in range(n)]
    X = pd.DataFrame(X, index=idx, columns=cols)
    y.index = batch.index = idx
    return quantile_normalize(batch_correct(X, batch, y, "zscore")), y, batch


# ----------------------------------------------------------------------------- modelling
def new_model(seed, **params):
    return XGBClassifier(eval_metric="logloss", random_state=seed, n_jobs=-1, **params)


def evaluate(model, X, y):
    p = model.predict_proba(X)[:, 1]
    yhat = (p >= 0.5).astype(int)
    m = dict(accuracy=accuracy_score(y, yhat), precision=precision_score(y, yhat, zero_division=0),
             recall=recall_score(y, yhat, zero_division=0), f1=f1_score(y, yhat, zero_division=0),
             auc=roc_auc_score(y, p))
    return m, confusion_matrix(y, yhat), p


def shap_importance(model, X):
    sv = shap.TreeExplainer(model).shap_values(X)
    if isinstance(sv, list):
        sv = sv[1]
    return sv, pd.Series(np.abs(sv).mean(axis=0), index=X.columns).sort_values(ascending=False)


def honest_cv(X, y, params, k, seed, repeats=3):
    """Gene selection by SHAP happens INSIDE each training fold -> no test leakage."""
    rows = []
    for tr, te in RepeatedStratifiedKFold(n_splits=5, n_repeats=repeats, random_state=seed).split(X, y):
        Xtr, Xte, ytr, yte = X.iloc[tr], X.iloc[te], y.iloc[tr], y.iloc[te]
        full = new_model(seed, **params).fit(Xtr, ytr)
        top = shap_importance(full, Xtr)[1].index[:k]
        small = new_model(seed, **params).fit(Xtr[top], ytr)
        rows.append(dict(acc_all=evaluate(full, Xte, yte)[0]["accuracy"], auc_all=evaluate(full, Xte, yte)[0]["auc"],
                         acc_topk=evaluate(small, Xte[top], yte)[0]["accuracy"],
                         auc_topk=evaluate(small, Xte[top], yte)[0]["auc"]))
    return pd.DataFrame(rows)


def diff_expr(X, y):
    a, b = X[y == 1], X[y == 0]
    _, p = stats.ttest_ind(a, b, equal_var=False, axis=0)
    de = pd.DataFrame({"mean_diff": a.mean() - b.mean(), "pval": p}, index=X.columns)
    de["FDR"] = multipletests(de["pval"].fillna(1), method="fdr_bh")[1]
    de["direction"] = np.where(de["mean_diff"] > 0, "up in positive class", "down in positive class")
    return de


def enrichment(genes, out):
    try:
        import gseapy as gp
        res = gp.enrichr(gene_list=list(genes), organism="human", outdir=None,
                         gene_sets=["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"]).results
        res = res[res["P-value"] < 0.05].sort_values("P-value")
        res.to_csv(out / "enrichment_enrichr.csv", index=False)
        print(f"Enrichr: {len(res)} terms with p<0.05 -> enrichment_enrichr.csv")
    except Exception as e:
        print(f"[skip] Enrichr failed (needs internet + gseapy): {e}")


# ----------------------------------------------------------------------------- plots
def make_plots(out, X, y, batch, Xtr, sv, model, Xte, yte, cm, k):
    pcs = PCA(n_components=2, random_state=0).fit_transform(X.fillna(0))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for name, grp, a in (("class", y, ax[0]), ("batch", batch, ax[1])):
        for v in sorted(grp.unique(), key=str):
            m = (grp == v).values
            a.scatter(pcs[m, 0], pcs[m, 1], s=14, label=str(v))
        a.set_title(f"PCA by {name}"); a.legend()
    fig.tight_layout(); fig.savefig(out / "pca_qc.png", dpi=150); plt.close(fig)

    for kind in ("bar", "dot"):
        shap.summary_plot(sv, Xtr, plot_type=kind, max_display=k, show=False)
        plt.tight_layout(); plt.savefig(out / f"shap_{kind}.png", dpi=150, bbox_inches="tight"); plt.close()

    p = model.predict_proba(Xte)[:, 1]
    fpr, tpr, _ = roc_curve(yte, p)
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    ax[0].plot(fpr, tpr, label=f"AUC={roc_auc_score(yte, p):.2f}"); ax[0].plot([0, 1], [0, 1], "k--")
    ax[0].set_xlabel("FPR"); ax[0].set_ylabel("TPR"); ax[0].set_title("ROC (test set)"); ax[0].legend()
    ax[1].imshow(cm, cmap="Blues")
    for (i, j), v in np.ndenumerate(cm):
        ax[1].text(j, i, int(v), ha="center", va="center")
    ax[1].set_xticks([0, 1]); ax[1].set_yticks([0, 1]); ax[1].set_title("Confusion matrix")
    ax[1].set_xlabel("predicted"); ax[1].set_ylabel("actual")
    fig.tight_layout(); fig.savefig(out / "roc_confusion.png", dpi=150); plt.close(fig)


# ----------------------------------------------------------------------------- main
def run(X, y, batch, out, k, seed, do_enrichr, fast=False, cv_repeats=3):
    out.mkdir(parents=True, exist_ok=True)
    print(f"\nFinal matrix: {X.shape[0]} samples x {X.shape[1]} genes | class counts {y.value_counts().to_dict()}")

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, stratify=y, random_state=seed)
    if fast:
        params = FAST_PARAMS.copy()
        print("Fast-mode params:", params)
        model = new_model(seed, **params).fit(Xtr, ytr)
    else:
        cv = min(3, int(ytr.value_counts().min()))
        gs = GridSearchCV(new_model(seed), GRID, cv=cv, scoring="roc_auc", n_jobs=1).fit(Xtr, ytr)
        params = gs.best_params_
        print("Best params:", params)
        model = gs.best_estimator_

    m_all, cm_all, _ = evaluate(model, Xte, yte)
    sv, imp = shap_importance(model, Xtr)               # SHAP from TRAIN data only
    top = list(imp.index[:k])
    small = new_model(seed, **params).fit(Xtr[top], ytr)
    m_top, cm_top, _ = evaluate(small, Xte[top], yte)

    summary = pd.DataFrame({"all_genes": m_all, f"top{k}_genes": m_top}).round(3)
    print("\nHold-out test metrics:\n", summary.to_string())
    print("Confusion matrix (all genes):\n", cm_all)

    hc = honest_cv(X, y, params, k, seed, repeats=cv_repeats)
    hc_sum = hc.agg(["mean", "std"]).round(3)
    print("\nLeakage-free repeated 5-fold CV (SHAP selection inside each fold):\n", hc_sum.to_string())

    de = diff_expr(X, y)
    table = pd.DataFrame({"mean_abs_SHAP": imp.iloc[:k]}).join(de)
    print(f"\nTop {k} genes by mean |SHAP|:\n", table.round(4).to_string())

    summary.to_csv(out / "metrics_holdout.csv"); hc.to_csv(out / "metrics_cv_folds.csv")
    table.to_csv(out / "top_genes_shap_de.csv"); imp.to_csv(out / "all_genes_shap_importance.csv")
    make_plots(out, X, y, batch, Xtr, sv, model, Xte, yte, cm_all, k)
    if do_enrichr:
        enrichment(top, out)
    print(f"\nAll outputs saved in: {out.resolve()}")
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", choices=PRESETS)
    ap.add_argument("--gse", nargs="+", help="custom GEO accession(s)")
    ap.add_argument("--label-col", help="only look for pos/neg words in metadata columns whose name contains this")
    ap.add_argument("--pos", nargs="+"); ap.add_argument("--neg", nargs="+"); ap.add_argument("--drop", nargs="*", default=None)
    ap.add_argument("--demo", action="store_true", help="synthetic data, tests the pipeline offline")
    ap.add_argument("--inspect", action="store_true", help="print GEO sample metadata and exit")
    ap.add_argument("--top-k", type=int, default=20)
    ap.add_argument("--batch-method", choices=["combat", "zscore", "none"], default="combat")
    ap.add_argument("--mitocarta", help="MitoCarta3.0 .xls or one-symbol-per-line text file")
    ap.add_argument("--cache", default="geo_cache"); ap.add_argument("--outdir")
    ap.add_argument("--seed", type=int, default=42); ap.add_argument("--no-enrichr", action="store_true")
    ap.add_argument("--fast", action="store_true", help="skip grid search and use fixed, reproducible XGBoost settings")
    ap.add_argument("--cv-repeats", type=int, default=3, help="repeats for stratified 5-fold CV (default: 3)")
    a = ap.parse_args()
    if a.cv_repeats < 1:
        ap.error("--cv-repeats must be at least 1")

    if a.demo:
        X, y, batch = demo_dataset(a.seed)
        name = "demo"
    else:
        cfg = dict(PRESETS.get(a.dataset, {}))
        gse = a.gse or cfg.get("gse")
        pos, neg = a.pos or cfg.get("pos"), a.neg or cfg.get("neg")
        drop = a.drop if a.drop is not None else cfg.get("drop", [])
        per = {} if (a.pos or a.neg) else cfg.get("per_gse", {})
        if not (gse and ((pos and neg) or per)):
            ap.error("give --dataset, or --gse with --pos and --neg (or use --demo)")
        Path(a.cache).mkdir(exist_ok=True)
        if a.inspect:
            for g in gse:
                load_gse(g, a.cache, inspect=True)
            print("\nExpected class counts:", cfg.get("expected", "see paper"))
            return
        print("Expected class counts:", cfg.get("expected", "see paper"))
        X, y, batch = build_dataset(gse, pos, neg, drop, a.cache, a.batch_method, a.mitocarta,
                                    a.label_col or cfg.get("label_col"), per)
        name = a.dataset or "_".join(gse)

    run(X, y, batch, Path(a.outdir or f"results_{name}"), a.top_k, a.seed,
        not (a.no_enrichr or a.demo), fast=a.fast, cv_repeats=a.cv_repeats)


if __name__ == "__main__":
    main()

