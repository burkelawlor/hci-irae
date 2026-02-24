# HCI irAE

Spatial transcriptomics analyses of immune-related adverse events (irAE) using Xenium data.

## About

Notebooks are layed out in the following order:

1. **01_preprocessing_qc.ipynb** — Load Xenium data, preprocessing, QC  
2. **02a_general_annotation_all_samples.ipynb** / **02b_general_annotation_subset.ipynb** — General cell-type annotation  
3. **03_immune_annotation.ipynb** — Immune cell annotation  
4. **04_lymphoid_regions.ipynb** — Lymphoid region analysis  

Shared utilities live in `utils/` (plotting, processing, cell types). Plots are written to `figures/`. Input data is expected under `data/` (not in repo).
