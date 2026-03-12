import math
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import squidpy as sq
import scanpy as sc

def spatial_plot_cell_types(
    adata,
    sample_name,
    ct_col,
    ncols=4,
    figsize=None,
    palette=None,
    size=10,
    save=False,
    na_label="NA",          # label to use for missing values plot
    other_label="Other",    # label for non-target cells in each panel
):
    adata_sample = adata[adata.obs["sample_name"] == sample_name].copy()

    # Identify cell types, including NA as its own category/panel
    ct_series = adata_sample.obs[ct_col]
    has_na = ct_series.isna().any()

    # Non-missing types (sorted for stable ordering)
    cts = list(np.sort(ct_series.dropna().unique().astype(str)))
    if has_na:
        cts.append(na_label)

    n_ct = len(cts)
    nrows = math.ceil(n_ct / ncols)

    if figsize is None:
        figsize = (4 * ncols, 4 * nrows)

    fig, ax = plt.subplots(nrows, ncols, figsize=figsize, tight_layout=True)
    ax = np.array(ax).reshape(-1)  # flatten safely for 1 row/col cases

    for i, ct in enumerate(cts):
        # Build a per-panel mask (special handling for NA)
        if ct == na_label:
            mask = ct_series.isna()
            panel_label = na_label
        else:
            # compare as strings to be robust to mixed dtypes
            mask = ct_series.astype(str) == ct
            panel_label = ct

        # Make an object/string series to avoid dtype promotion issues
        is_ct = np.where(mask.to_numpy(), panel_label, other_label).astype(object)
        adata_sample.obs["is_ct"] = is_ct

        # Background: all cells (no color)
        sq.pl.spatial_scatter(
            adata_sample,
            library_id="spatial",
            shape=None,
            ax=ax[i],
            na_color="lightgray",
            size=size,
        )

        # Foreground: highlight the current ct (colored by is_ct)
        adata_fg = adata_sample[mask].copy()
        if adata_fg.n_obs > 0:
            if palette is not None and panel_label in palette:
                sq.pl.spatial_scatter(
                    adata_fg,
                    library_id="spatial",
                    shape=None,
                    color="is_ct",
                    ax=ax[i],
                    size=size,
                    palette=ListedColormap([palette[panel_label]]),
                )
            else:
                sq.pl.spatial_scatter(
                    adata_fg,
                    library_id="spatial",
                    shape=None,
                    color="is_ct",
                    ax=ax[i],
                    size=size,
                )

        ax[i].set_title(f"{panel_label}")
        ax[i].invert_yaxis()

        leg = ax[i].get_legend()
        if leg is not None:
            leg.remove()
        ax[i].set_axis_off()

    # Delete unused axes
    for j in range(n_ct, len(ax)):
        fig.delaxes(ax[j])

    plt.suptitle(f"{sample_name}")
    plt.show()

    if save:
        fig.savefig(f"figures/spatial_plots/{sample_name}_{ct_col}_individual.png", bbox_inches="tight")




def feature_plots_from_marker_genes(adata, marker_genes_dict, save=False, prefix=None):
    
    for ct in marker_genes_dict:
        print(f'{ct.upper()}:')
        
        ax = sc.pl.umap(adata, color=marker_genes_dict[ct], wspace=0.1, show=False)
        plt.suptitle(ct)
        
        if save:
            if isinstance(ax, list):
                ax[0].figure.savefig(f"./figures/feature_plots/{prefix}_umap_{ct}.png", bbox_inches="tight")
            else:
                ax.figure.savefig(f"./figures/feature_plots/{prefix}_umap_{ct}.png", bbox_inches="tight")

        plt.show()
    