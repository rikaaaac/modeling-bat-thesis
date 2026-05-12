import os
import scanpy as sc
import anndata as ad
import numpy as np

# Set the path to the folder containing .h5 files
h5_path = '/Users/rikac/Downloads/h5/'

# Initialize a dictionary to hold the data
bc_matrix_dict = {}

# Iterate through files in the directory
for filename in os.listdir(h5_path):
    if filename.startswith('JJ0') and filename.endswith('.h5'):
        samplename = filename.replace('_filtered_feature_bc_matrix.h5', '')

        # Load data from 10x-formatted h5 file
        adata = sc.read_10x_h5(h5_path + filename, gex_only=True)
        adata.var_names_make_unique()
        adata.obs['sample'] = samplename

        bc_matrix_dict[samplename] = adata

# Combine selected samples into one AnnData object
samples_to_combine = ["JJ003", "JJ004", "JJ005"]
adatas_to_concat = [bc_matrix_dict[sample] for sample in samples_to_combine if sample in bc_matrix_dict]

if len(adatas_to_concat) == 0:
    raise ValueError("No valid samples found to combine.")

our_adata = ad.concat(adatas_to_concat, join='outer')
our_adata.obs_names_make_unique()

# Show number of cells per sample
print(our_adata.obs['sample'].value_counts())

# # filter to keep only 15,000 cells and 2000 genes

# # Sample cells (rows)
# if our_adata.n_obs > 15000:
#     our_adata = our_adata[np.random.choice(our_adata.obs_names, size=15000, replace=False), :]

# # Sample genes (columns)
# if our_adata.n_vars > 2000:
#     our_adata = our_adata[:, np.random.choice(our_adata.var_names, size=2000, replace=False)]

# our_adata.write_h5ad('/Users/rikac/Documents/test/sample_adata.h5ad')

