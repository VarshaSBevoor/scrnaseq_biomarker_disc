# Project 1: scRNA-seq Biomarker Discovery with ML

## What this project is
A single-cell RNA-seq analysis pipeline that goes beyond differential gene expression to use ML + SHAP for cell-type-aware biomarker discovery. Built for two purposes: (1) GitHub portfolio visibility, (2) deep personal learning so I can speak confidently about every decision in interviews.

## Dataset
- **Kang et al. 2018 Lupus PBMC dataset** (stimulated vs control PBMCs)
- Available via CELLxGENE or scanpy datasets
- ~14,000 cells, 8 donors, IFN-β stimulated vs control condition

## Tech stack
- `scanpy` — preprocessing, clustering, annotation
- `scikit-learn` — ML classifier (start with logistic regression, then try random forest)
- `SHAP` — model interpretability, biomarker candidate extraction
- `matplotlib`, `seaborn` — visualization
- Python only

## Project phases
1. **QC + preprocessing** — doublet detection, filtering, normalization, HVG selection
2. **Clustering + cell type annotation** — Leiden clustering, marker-based annotation
3. **ML classifier** — disease (stimulated) vs control, trained per cell type
4. **SHAP analysis** — extract top predictive genes per cell type → biomarker candidates
5. **Biological validation** — cross-reference candidates with known IFN-response and lupus literature

## Learning goals (not just implementation)
- Understand WHY each preprocessing step exists, not just how to run it
- Understand what SHAP values mean mathematically and biologically
- Be able to explain every parameter choice (e.g. n_neighbors, resolution in Leiden)
- Know how to talk about limitations (class imbalance, batch effects, overfitting risk)

## Interview narrative
"I went beyond standard DGE — I trained a cell-type-specific classifier to distinguish stimulated from control cells, then used SHAP to identify which genes drove those predictions. That gives you interpretable, ML-derived biomarker candidates that you can validate against known biology."

## Working style
- Explain concepts as we build — understanding over speed
- Flag when I should stop and read something before proceeding
- Suggest meaningful commit messages that document learning, not just working code
- Commit typed prefixes: feat, fix, explore, data, docs, refactor

## Current status
- [ ] Project scaffolded
- [ ] Dataset downloaded
- [ ] QC complete
- [ ] Clustering + annotation complete
- [ ] ML classifier trained
- [ ] SHAP analysis done
- [ ] Biological validation done
- [ ] README written
