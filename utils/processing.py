import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad


def propogate_subset_labels(parent_adata, subset_adata, key_added, key_to_add, merge_on):
    """
    Propogate labels from subset_adata.obs[key_to_add] to parent_adata.obs[key_added].

    Args:
        parent_adata: AnnData object to add the labels to
        subset_adata: AnnData object to get the labels from
        key_added: Name of the new column to add to parent_adata.obs
        key_to_add: Name of the column in subset_adata.obs to add to parent_adata.obs
        merge_on: Name of the column in subset_adata.obs and parent_adata.obs to merge on
    """

    adata_new = parent_adata.copy()

    if key_added not in adata_new.obs.columns:
        adata_new.obs[key_added] = np.nan

    adata_new.obs = adata_new.obs.merge(subset_adata.obs[[merge_on, key_to_add]], how="left", on=merge_on)
    adata_new.obs[[key_added, key_to_add]] = adata_new.obs[[key_added, key_to_add]].astype('object')
    
    mask = adata_new.obs[key_to_add].notna()    
    adata_new.obs.loc[mask, key_added] = adata_new.obs.loc[mask, key_to_add]
    adata_new.obs[key_added] = adata_new.obs[key_added].astype("category")
    adata_new.obs.drop(columns=[key_to_add], inplace=True)
    
    return adata_new