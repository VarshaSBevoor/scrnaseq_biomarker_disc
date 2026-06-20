# scRNA-seq Biomarker Discovery with ML + SHAP

Single-cell RNA-seq analysis pipeline that goes beyond standard differential gene expression to identify **interpretable, cell-type-specific biomarker candidates** using machine learning and SHAP.

## Dataset

**Kang et al. 2018** — IFN-β stimulated vs control PBMCs (GSE96583)
- ~28,500 cells, 8 donors, 8 cell types
- Condition: interferon-beta stimulated vs untreated control

## Approach

Standard differential gene expression tests genes independently and ignores cell type context. This pipeline instead:

1. Trains a **per-cell-type logistic regression classifier** to predict stimulation status
2. Applies **SHAP** to explain which genes drove each cell's prediction
3. Identifies **consensus biomarker candidates** — genes that are top predictors across multiple cell types independently

This gives predictive, interpretable, cell-type-specific biomarkers rather than simple statistical associations.

## Pipeline

| Phase | Notebook | Description |
|-------|----------|-------------|
| 1 | `01_qc_preprocessing` | QC filtering, normalization, HVG selection |
| 2 | `02_clustering` | PCA, Leiden clustering, UMAP visualization |
| 3 | `03_ml_classifier` | Per-cell-type logistic regression (donor-aware split) |
| 4 | `04_shap_analysis` | SHAP feature importance, biomarker extraction |
| 5 | `05_biological_validation` | Cross-reference with IFN-β literature |

## Results

Classifiers achieved **F1 > 0.95** for 7/8 cell types on held-out donors.

Top consensus biomarker candidates (appearing as top SHAP features across multiple cell types):

| Gene | Cell types | Known role |
|------|-----------|------------|
| ISG15 | 8/8 | Hallmark IFN-stimulated gene; proposed lupus biomarker |
| IFIT3 | 7/8 | Interferon-induced; blocks viral replication |
| LY6E | 7/8 | IFN-induced; restricts viral entry |
| ISG20 | 7/8 | IFN-stimulated exonuclease; degrades viral RNA |
| IFIT1 | 6/8 | Interferon-induced; inhibits viral translation |
| IRF7 | 4/8 | Master transcription factor driving IFN production |

The ML model, trained with no biological knowledge, independently rediscovered the known IFN-β response signature — providing strong biological validation of the approach.

## Tech Stack

- `scanpy` — preprocessing, clustering, visualization
- `scikit-learn` — logistic regression classifier
- `shap` — model interpretability
- `matplotlib`, `seaborn` — visualization

## Setup

```bash
conda create -n scrnaseq-biomarker python=3.11
conda activate scrnaseq-biomarker
pip install scanpy scikit-learn shap scrublet jupyter matplotlib seaborn
```

Data is downloaded automatically from GEO (GSE96583) when running notebook 01.
