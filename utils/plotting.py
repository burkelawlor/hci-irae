import math
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import squidpy as sq
import scanpy as sc

def spatial_plot_cell_types(adata, sample_name, ct_col, ncols=4, figsize=None, palette=None, size=10, save=False):
    
    adata_sample = adata[adata.obs['sample_name'] == sample_name].copy()
    n_ct = len(adata_sample.obs[ct_col].cat.categories)
    nrows = math.ceil(n_ct / ncols)

    fig, ax = plt.subplots(nrows, ncols, figsize=figsize, tight_layout=True)
    ax = ax.flatten()
    nrows = math.ceil(n_ct / ncols)
    for i, ct in enumerate(adata_sample.obs[ct_col].cat.categories):

        adata_sample.obs['is_ct'] = np.where(adata_sample.obs[ct_col] == ct, ct, 'Other')

        sq.pl.spatial_scatter(adata_sample, library_id="spatial", shape=None, ax=ax[i], na_color='lightgray', size=10)
        
        if palette is not None:
            sq.pl.spatial_scatter(adata_sample[adata_sample.obs['is_ct'] == ct], library_id="spatial", shape=None, color='is_ct', ax=ax[i], size=size, palette=ListedColormap([palette[ct]]))
        else:
            sq.pl.spatial_scatter(adata_sample[adata_sample.obs['is_ct'] == ct], library_id="spatial", shape=None, color='is_ct', ax=ax[i], size=size)
            
        ax[i].set_title(f'{ct}')
        ax[i].invert_yaxis()
        ax[i].get_legend().remove()
        ax[i].set_axis_off()

    # Delete unused axes
    for j in range(n_ct, len(ax)):
        fig.delaxes(ax[j])

    plt.suptitle(f'{sample_name}')
    plt.show()

    if save:
        ax[0].figure.savefig(f'figures/spatial_plots/{sample_name}_{ct_col}_individual.png', bbox_inches='tight')




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
    