import pandas as pd
import scanpy as sc
import anndata as ad

def propogate_subset_labels(parent_adata, subset_adata, key_added, key_to_add, merge_on):
    """
    Propogate labels from subset_adata.obs[key_to_add] to parent_adata.obs[key_added].
    """
    adata_new = parent_adata.copy()
    adata_new.obs = adata_new.obs.merge(subset_adata.obs[[merge_on, key_to_add]], how="left", on=merge_on)
    adata_new.obs[[key_added, key_to_add]] = adata_new.obs[[key_added, key_to_add]].astype('object')
    mask = adata_new.obs[key_to_add].notna()    
    adata_new.obs.loc[mask, key_added] = adata_new.obs.loc[mask, key_to_add]
    adata_new.obs[key_added] = adata_new.obs[key_added].astype("category")
    adata_new.obs.drop(columns=[key_to_add], inplace=True)
    return adata_new