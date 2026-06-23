# Learning Log — scRNA-seq Biomarker Discovery

---

## Session 1 — Phase 1: QC + Preprocessing (setup & data exploration)

### Q: Why are the count matrix and metadata stored as separate files?

Sparsity is part of it (scRNA-seq matrices are ~90% zeros), but the deeper reason is that they come from **different steps and different tools**:
- The **count matrix** comes from the sequencer + Cell Ranger (10x's alignment pipeline). It only knows "cell barcode X had Y reads mapping to gene Z." It knows nothing about biology.
- The **metadata** (condition, donor, cell type) comes from downstream analysis — demultiplexing to assign donor identity, and cell typing to label each cell. This is human interpretation layered on top.

So they're separate because they're produced at different stages, potentially by different people.

---

### Q: Where does the data come from?

**NCBI GEO** (Gene Expression Omnibus) — the standard public repository for genomics datasets. Every published scRNA-seq paper deposits its data here with an accession number.

- Kang 2018 accession: **GSE96583**
- URL: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583`

**Files we use (batch 2 = the PBMC IFN-β experiment):**
| File | What it contains |
|------|-----------------|
| `GSE96583_RAW.tar` | Raw count matrices (MTX format, one per condition) |
| `GSE96583_batch2.genes.tsv.gz` | Gene names (Ensembl ID + symbol) |
| `GSE96583_batch2.total.tsne.df.tsv.gz` | Cell metadata: condition, donor, cell type, doublet flag |

---

### Q: What is AnnData?

Scanpy's core data structure. It keeps everything linked:
- `.X` — the count matrix (cells × genes)
- `.obs` — cell-level annotations (rows), indexed by cell barcode
- `.var` — gene-level annotations (columns), indexed by gene name
- `.layers` — alternative matrices (e.g., raw counts alongside normalized counts)

The key insight: rows in `.obs` are **cells**, not genes. Each row is one cell, indexed by its barcode (e.g., `AAACATACAATGCC-1`).

---

### Dataset snapshot (from `meta.head()` and `value_counts()`):

**Shape:** 29,065 cells × 7 metadata columns

**Condition balance:**
```
ctrl    14,619
stim    14,446
```
Nearly 50/50 — ideal for a classifier. Severe imbalance (e.g., 90/10) would require techniques like class weighting or oversampling.

**Cell types:**
```
CD4 T cells          12,033   (~41%)
CD14+ Monocytes       6,447
B cells               2,880
CD8 T cells           2,634
NK cells              2,330
FCGR3A+ Monocytes     1,914
Dendritic cells          472
Megakaryocytes           346   (~1%)
```
CD4 T cells dominate; Megakaryocytes are rare. This matters in Phase 3 — per-cell-type classifiers will have much less training data for rare cell types.

**Donors:**
```
1015    5,841
1488    5,290
1256    4,726
1016    4,178
1244    4,005
101     2,390
1039    1,349
107     1,286
```
8 donors total. Donor IDs matter because of **batch effects** — a model might learn donor identity instead of biological signal (stimulated vs control). Fix: use a donor-aware train/test split in Phase 3 (hold out entire donors, not random cells).

---

---

### Q: What is MTX format?

MTX (Matrix Market Exchange) is a format for **sparse matrices**. Instead of storing a full grid (most of which would be zeros), it only stores non-zero values as `(row, col, value)` triplets:

```
% comment
35635 14619 3421847    ← (n_genes, n_cells, n_nonzero_entries)
1 1 4                  ← gene 1, cell 1 has count 4
1 5 2                  ← gene 1, cell 5 has count 2
```

Three files always go together:
- `.mtx.gz` — the sparse counts
- `barcodes.tsv.gz` — cell identifiers (one per column)
- `genes.tsv.gz` — gene names (one per row)

**Important:** MTX stores genes × cells, but AnnData expects cells × genes — so you always transpose (`.T`) when loading.

---

---

### Q: Why use Ensembl IDs instead of gene symbols as the index?

Gene symbols are human-readable but not unique — two different genomic loci can share the same symbol, and symbols change over time as annotation databases update. Ensembl IDs (`ENSG00000...`) are stable and always unique.

**Rule:** Use Ensembl IDs as the var index, keep symbols as a column for readability.

---

### Q: Why save raw counts in `adata.layers["counts"]` before normalizing?

Once you normalize and log-transform `adata.X`, the original counts are gone. You need them later for:
1. **HVG selection** — the variance modeling algorithm expects raw counts
2. **SHAP analysis** — raw counts let you cross-reference interpretations

`adata.layers` stores multiple versions of the matrix simultaneously. `adata.layers["counts"]` = original, `adata.X` = normalized.

---

### Where we stopped:
Just saved raw counts to `adata.layers["counts"]`. AnnData object is fully assembled:
- 29,065 cells × 35,635 genes
- `.obs` columns: `sample`, `donor`, `condition`, `cell_type`, `multiplets`
- `.layers`: `counts`

**Next step:** Visualizing QC distributions and choosing filter thresholds.

---

### Q: What makes a cell "low quality" in scRNA-seq?

Three things to check, each with a metric:

| Problem | Metric | Direction |
|---------|--------|-----------|
| Empty droplet (ambient RNA, no real cell) | `n_genes_by_counts`, `total_counts` | Too **low** |
| Doublet (two cells in one droplet) | `n_genes_by_counts`, `total_counts` | Too **high** |
| Dying/damaged cell | `pct_counts_mt` | Too **high** |

**Why high mito % = dying cell:** When a cell is stressed, the cytoplasm (containing nuclear mRNA) leaks out but mitochondria are retained inside the membrane. You end up sequencing mostly mitochondrial transcripts.

**Why doublets have high counts:** Two cells pooled = ~2x the genes and UMIs of a normal single cell.

---

### Q: Why did `gene_symbol` disappear from `adata_combined.var` after `ad.concat`?

`ad.concat` only keeps var columns that are consistent across all objects being merged. If there's any issue with alignment it drops them. Fix: re-assign after concat:
```python
adata_combined.var["gene_symbol"] = genes["gene_symbol"].values
```

---

### Code so far (session 3):
```python
adata_combined.obs_names_make_unique()  # fix 313 duplicate barcodes from two 10x runs
adata_combined.obs = adata_combined.obs.join(meta[["ind", "stim", "cell", "multiplets"]])
adata_combined.obs.rename(columns={"ind": "donor", "stim": "condition", "cell": "cell_type"}, inplace=True)
adata_combined.layers["counts"] = adata_combined.X.copy()
adata_combined.var["gene_symbol"] = genes["gene_symbol"].values
adata_combined.var["mt"] = adata_combined.var["gene_symbol"].str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_combined, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
```

**13 mitochondrial genes** detected (expected for human — mitochondrial genome encodes exactly 13 proteins).

---

---

### Q: Why normalize scRNA-seq data?

Each cell is sequenced to a different depth — one cell might have 5,000 total UMIs, another 50,000. Without normalization, highly-sequenced cells would appear to express every gene more, purely as a technical artifact.

**Two-step normalization:**
1. `normalize_total(target_sum=1e4)` — scales each cell so its total counts sum to 10,000 (counts per 10k)
2. `log1p` — compresses the dynamic range. The `+1` keeps zeros as zeros. Converts multiplicative differences → additive, which most statistical methods expect.

**Example:** Gene with 151 raw counts in a cell with 9,677 total:
- After normalize_total: `151 / 9677 × 10000 = 156`
- After log1p: `log(156 + 1) = 5.06`

---

### Q: What are Highly Variable Genes (HVGs) and why select them?

Housekeeping genes (e.g., GAPDH) are expressed similarly in every cell — not informative. HVGs vary meaningfully across cells, capturing differences between cell types and states. We select the top 2,000.

**Key parameters in `sc.pp.highly_variable_genes`:**
- `flavor="seurat_v3"` — fits a negative binomial model to estimate biological vs technical variance
- `layer="counts"` — uses raw counts because the model is designed for integer counts, not normalized values
- `batch_key="donor"` — computes HVGs within each donor separately, avoiding genes that look variable just because of donor-to-donor differences
- `subset=False` — keeps all genes in the object but flags the top 2,000 in `var["highly_variable"]`

---

### Q: Why set a minimum gene filter of 200 (not 0)?

Cells with very few detected genes are likely empty droplets that captured ambient RNA, not a real cell. 200 is a conservative threshold — the real cell distribution in this dataset starts around 500 genes.

---

### Q: Why filter genes expressed in fewer than 10 cells?

Statistical reliability — a gene detected in fewer than 10 cells gives almost no signal. Any model or test would be unreliable. Also practical: sparse genes bloat the matrix with zeros.

---

## Phase 1 complete ✓

**Final object:** 28,557 cells × 13,727 genes  
**Saved to:** `data/raw/kang18_preprocessed.h5ad`

**Filters applied:**
- Cells: 200 ≤ n_genes ≤ 2500, total_counts ≤ 10,000, no missing cell_type
- Genes: expressed in ≥ 10 cells
- Note: mito filter skipped — GEO data had no mitochondrial counts (filtered during alignment)

---

## Phase 2: Clustering + Cell Type Annotation

### Pipeline: PCA → neighbors → Leiden → UMAP

**Why PCA first?**
13,727 genes = 13,727 dimensions. Most algorithms break at that scale. PCA compresses to ~50 dimensions while keeping biological variation. Only run on the 2,000 HVGs — noise genes would dilute the signal.

**Why scale before PCA?**
PCA finds directions of maximum variance. A gene with range 0-5000 would dominate over one with range 0-0.5. Scaling (mean=0, std=1) puts every gene on equal footing. `max_value=10` clips extreme outliers.

**What is a principal component?**
An axis of variation capturing genes that move together. PC1 = biggest source of variation (e.g. T cells vs monocytes). PC2 = next biggest. Look at the variance ratio plot elbow to choose how many PCs to use. In this dataset: elbow at ~PC15-20, so we use 20.

**What is the neighbor graph?**
After PCA, each cell is a point in 20D space. `sc.pp.neighbors(n_pcs=20, n_neighbors=15)` connects each cell to its 15 most similar cells. Leiden clustering then finds densely connected groups. Without the graph, the clustering algorithm has no notion of similarity.

**Leiden clustering (`resolution=0.5`):**
- Higher resolution → more, finer clusters
- Lower resolution → fewer, coarser clusters
- Result: **16 clusters**

**UMAP:**
- Reduces to 2D for visualization only (not used for computation)
- Preserves local structure — similar cells land close together
- Caveat: distances between clusters are not quantitatively meaningful

### Key observations from UMAPs:
- Leiden clusters ≈ paper's cell type labels → clustering recovered real biology without supervision
- Some cell types split into 2 clusters → sub-populations (e.g. naive vs memory T cells) — biologically meaningful
- **Ctrl and stim cells are distinctly separated** → IFN-β stimulation causes a large transcriptional shift, clearly visible in UMAP space. This is why the ML classifier will work.
- Each cluster is dominated by one cell type → clean, interpretable clustering

### Why train the classifier per cell type (Phase 3)?
Different cell types respond differently to IFN-β. The genes distinguishing ctrl from stim in T cells are different from those in monocytes. Training per cell type gives you cell-type-specific biomarker candidates.

---

## Phase 3: ML Classifier (stimulated vs control, per cell type)

### Why per cell type?
Training one classifier on all cells would confound cell type identity with stimulation response — it would learn "T cell in ctrl" vs "monocyte in stim" rather than the actual biological signal. Training per cell type asks: *within this cell type, what genes distinguish ctrl from stim?*

### Input features
- Used 2,000 HVGs (not all 13,727 genes) — noise genes dilute signal and cause overfitting
- `X = np.array(subset.X)` — log-normalized expression values
- `y` — binary label: stim=1, ctrl=0

### Why donor-aware train/test split?
Random splitting lets the same donor's cells appear in both train and test. The model could learn donor-specific quirks and still get high accuracy — that's **data leakage**. Holding out entire donors tests whether the model generalizes to *new* donors, which is the real question.

### Model: Logistic Regression
- `C=1.0` — regularization parameter (inverse of strength). Low C = strong regularization = simpler model. With 2,000 features and relatively few cells, regularization prevents overfitting.
- `max_iter=1000` — gives the optimizer enough iterations to converge

### Understanding model evaluation metrics

**Confusion matrix concepts:**
- **True Positive (TP):** model predicted stim, cell is stim ✓
- **True Negative (TN):** model predicted ctrl, cell is ctrl ✓
- **False Positive (FP):** model predicted stim, cell is actually ctrl ✗
- **False Negative (FN):** model predicted ctrl, cell is actually stim ✗

**Precision** = TP / (TP + FP) → of all cells predicted as stim, how many actually are?
**Recall** = TP / (TP + FN) → of all actual stim cells, how many did we catch?
**F1** = harmonic mean of precision and recall = 2 × (P × R) / (P + R)

**Why F1 and not just accuracy?**
Accuracy can be misleading with imbalanced classes. If 90% of cells are ctrl, a model that always predicts ctrl gets 90% accuracy but is useless. F1 balances precision and recall, giving a fairer picture.

**Macro avg F1** = average F1 across both classes (ctrl and stim), treating them equally regardless of size.

### Results per cell type:
| Cell type | F1 score |
|-----------|----------|
| CD14+ Monocytes | 0.998 |
| FCGR3A+ Monocytes | 0.997 |
| B cells | 0.988 |
| CD8 T cells | 0.980 |
| NK cells | 0.970 |
| CD4 T cells | 0.963 |
| Dendritic cells | 0.959 |
| Megakaryocytes | 0.717 |

**Why Megakaryocytes are lower:** Only 346 cells total — after donor-aware split, train set is tiny. Not enough examples to learn a reliable boundary. Real limitation to flag.

**Why monocytes score highest:** Monocytes are primary responders to interferon signaling. IFN-β causes a strong, consistent transcriptional shift in monocytes across all donors.

**Why 100% (or near) isn't suspicious here:** IFN-β stimulation causes a massive, well-documented transcriptional response. The UMAP showed ctrl and stim completely separated — the biology is genuinely that strong. Verified no data leakage by confirming every test donor had both ctrl and stim cells.

---

## Phase 4: SHAP Analysis → Biomarker Candidates

### What is SHAP?
SHAP (SHapley Additive exPlanations) explains *why* a model made each prediction. For every cell, every gene gets a SHAP value:
- **Positive** → gene pushed prediction toward stim
- **Negative** → gene pushed prediction toward ctrl
- **Near zero** → gene didn't matter for this cell

Averaging absolute SHAP values across all cells gives a gene-level importance ranking — the top genes are your biomarker candidates.

**Advantage over model coefficients:** SHAP accounts for feature interactions and gives cell-level explanations, not just a global ranking.

### Reading a SHAP summary plot:
- **Red dots, positive SHAP** → gene is highly expressed in stim cells, drives stim prediction ✓ (biologically expected)
- **Blue dots, negative SHAP** → gene is lowly expressed in stim cells, drives ctrl prediction
- Wide spread of dots → gene is important across many cells
- Tight cluster near zero → gene doesn't matter much

### Top genes for CD14+ Monocytes:
CCL8, IFITM3, ISG15, APOBEC3A, IFIT2, LY6E, RSAD2, TNFSF10, IFIT3, ISG20, IFIT1, IDO1, CXCL10, CXCL11

All are known IFN-β response genes — the model rediscovered known biology with no prior biological knowledge.

### Consensus biomarker candidates (appearing across multiple cell types):
| Gene | # cell types |
|------|-------------|
| ISG15 | 8 (all cell types) |
| IFIT3 | 7 |
| LY6E | 7 |
| ISG20 | 7 |
| IFIT1 | 6 |
| IFIT2 | 4 |
| IRF7 | 4 |
| OAS1 | 3 |

**Key insight:** These genes (ISG = Interferon Stimulated Gene, IFIT = Interferon-Induced with Tetratricopeptide repeats, IRF = Interferon Regulatory Factor) are literally named after interferon. The ML model had no biological knowledge — it independently rediscovered the known IFN-β response signature.

**ISG15 in all 8 cell types** is the single strongest biomarker candidate — a universal IFN-β response gene across all PBMC populations.

### Main finding (interview narrative):
*"Unsupervised SHAP analysis of cell-type-specific classifiers converged on the known IFN-β response gene signature — ISG15, IFIT family, IRF7 — with ISG15 appearing as a top predictor in all 8 cell types. This validates that the models learned real biology rather than technical artifacts."*

### Saved outputs:
- `results/consensus_biomarkers.csv` — genes ranked by number of cell types they appear in
- `results/top_genes_<cell_type>.csv` — top 10 SHAP genes per cell type

---

## Phase 5: Biological Validation

### What is ISG15?
One of the most well-studied interferon-stimulated genes:
- Ubiquitin-like protein, massively upregulated within hours of IFN stimulation
- Tags other proteins (ISGylation) to modulate antiviral response
- Hallmark of interferon signaling — if ISG15 is up, IFN pathway is active
- Proposed biomarker in lupus — directly relevant to the Kang 2018 disease context

### Validation result:
10/20 top SHAP candidates have documented roles in IFN signaling. The top 8 by consensus (appearing in 6-8 cell types) are ALL validated known IFN-response genes. Genes appearing in only 1-2 cell types are more speculative.

---

## Why this approach is better than standard Differential Gene Expression (DGE)

**Standard DGE** asks: which genes are statistically different between ctrl and stim? Tests each gene independently, gives p-values and fold changes.

**This project adds three things DGE can't give you:**

1. **Cell-type specificity** — DGE on all cells mixes signals from different populations. Per-cell-type classifiers tell you which genes matter *within* each cell type separately.

2. **Predictive validation** — the genes don't just correlate with stimulation, they *predict* it with 96-99% F1. Prediction is a stronger claim than statistical association.

3. **Interpretability with SHAP** — you know not just which genes differ, but how much each gene contributed to each individual cell's classification. Cell-level resolution DGE can't provide.

**One-sentence interview answer:** *"DGE finds associations. My approach finds predictive, interpretable, cell-type-specific biomarkers — a stronger claim for translational relevance."*

---

## Project complete ✓

---

## Limitations (important for interviews)

Interviewers ask about limitations to test scientific maturity. Know these five:

**1. No independent validation dataset**
Held-out donors came from the same study. Ideally you'd test on a completely separate lupus dataset to show biomarkers generalize beyond Kang 2018.

**2. Megakaryocytes are underpowered**
Only 346 cells → F1 of 0.717. Results for rare cell types are unreliable. Need more cells to trust those biomarker candidates.

**3. Logistic regression assumes linearity**
Model assumes gene effects are additive and linear. Gene regulation is often non-linear (combinatorial effects, thresholds). Random forest or neural network might capture more complex patterns — though logistic regression's simplicity makes SHAP more interpretable.

**4. Batch effects not fully corrected**
Used `batch_key=donor` in HVG selection but no explicit batch correction (Harmony, scVI). Some signal might be donor-specific rather than biological.

**5. Correlation not causation**
SHAP identifies predictive genes — not causal drivers. ISG15 being upregulated doesn't mean it *causes* the IFN response; it may be a downstream marker. Functional validation (knockdown experiments) would be needed to establish causality.

---

## Interview Q&A

**Q: Walk me through your pipeline.**
Raw 10x data from GEO → QC filtering (empty droplets, doublets) → normalization → HVG selection → PCA + Leiden clustering → per-cell-type logistic regression classifier with donor-aware split → SHAP analysis → consensus biomarker candidates validated against IFN-β literature.

**Q: Why per cell type?**
Training on all cells confounds cell type identity with stimulation response. Per-cell-type classifiers ask: within this cell type, what genes distinguish ctrl from stim? That gives cell-type-specific biomarker candidates.

**Q: What is SHAP?**
Explains why a model made each prediction. For every cell, every gene gets a SHAP value — positive means it pushed toward stim, negative toward ctrl. Averaging absolute SHAP values across cells gives a gene importance ranking. The advantage over coefficients: cell-level resolution and accounts for feature interactions.

**Q: 99% F1 on monocytes — should I be suspicious?**
No — verified no data leakage (every test donor had both ctrl and stim cells). The UMAP showed ctrl and stim completely separated for monocytes. IFN-β has a massive, well-documented transcriptional effect on monocytes, which are primary interferon responders. The biology is genuinely that strong.

**Q: How is this better than DGE?**
DGE finds associations. This approach finds predictive, interpretable, cell-type-specific biomarkers. Three advantages: (1) cell-type specificity — DGE mixes signals across populations, (2) predictive validation — genes predict stimulation at 96-99% F1, not just correlate with it, (3) cell-level resolution via SHAP — DGE can't tell you how much each gene contributed to each individual cell's classification.

**Q: What are the limitations?**
No independent validation dataset, Megakaryocytes underpowered, logistic regression assumes linearity, batch effects not fully corrected, correlation not causation.

---

## TODO:
- [x] Phase 1: QC + preprocessing
- [x] Phase 2: Clustering + cell type annotation
- [x] Phase 3: ML classifier (stimulated vs control, per cell type)
- [x] Phase 4: SHAP analysis → biomarker candidates
- [x] Phase 5: Biological validation
- [x] README
- [ ] Practice explaining out loud
