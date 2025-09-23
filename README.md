# CellCensusNMF

Reproducible code & outputs for the manuscript figures.

**👉 See the rendered code + comments + images:**  
[Figures_code.md](Figures_code.md)

---

## What’s here

- **`Figures_code.Rmd`** – the source R Markdown with all code & inline commentary.  
- **`Figures_code.md`** – the knitted, GitHub-friendly output (shows code + PNGs).  
- **`fig/`** – PNGs used by the markdown (e.g., `Figure2.png`).  
---

## Quick start 

If already have the precomputed inputs, knit in **QUICK** mode to reproduce the figures without re-running long analyses.

## Full compute 

Run FULL mode to recompute results from raw/intermediate inputs (e.g., .h5ad, embeddings). This needs more dependencies/data.

rmarkdown::render("Figures_code.Rmd", params = list(mode = "full"))

- Required inputs for FULL mode (see [cell_census_datasets.py](cell_census_datasets.py)):

adata_obsm.csv

adata_obs.csv

adata_var_metadata.csv

adata_var_embeddings.csv

transfer_learning_data.h5ad

