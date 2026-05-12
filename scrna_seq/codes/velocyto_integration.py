# RNA velocity
import scanpy as sc
import numpy as np
import velocyto as vcy
import loompy as lp
import scipy.sparse as sp
import scvelo as scv
import cellrank as cr
import re
import os
os.environ["CELLRANK_DISABLE_JIT"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMBA_DISABLE_JIT"] = "1"

import multiprocessing as mp
mp.set_start_method("fork", force=True)

jj003_loom = '/Users/rikac/Documents/mansfield_lab/spring_25/scrna_seq_analysis/JJ003.loom'
jj004_loom = '/Users/rikac/Documents/mansfield_lab/spring_25/scrna_seq_analysis/JJ004.loom'
jj005_loom = '/Users/rikac/Documents/mansfield_lab/spring_25/scrna_seq_analysis/JJ005.loom'
file1 = '/Users/rikac/Documents/mansfield_lab/spring_25/scrna_seq_analysis/bc3306_final_adata.h5ad'



#######################################################################################
# INTEGRATING LOOM FILE INTO H5AD 

vlm = vcy.VelocytoLoom(jj003_loom)

with lp.connect(jj003_loom) as ds:
    cell_names = ds.ca["CellID"]
    print(cell_names[:10])

# normalize
vlm.normalize("S", size=True, log=True)

# plot
vlm.plot_fractions()

# test - integrat JJ003 loom into h5ad file
file = '/Users/rikac/Documents/project_lab/spring_25/scrna_seq_analysis/bc3306_mt_filter_sc_updated_2.h5ad'
adata = sc.read(file)

# make sure both loom and h5ad have the same barcode
loom_barcodes = [c.split(":")[1].rstrip("x") for c in cell_names]

subset = adata[adata.obs["sample"] == "JJ003"].copy() # 9124 × 18973
subset.obs["barcode"] = subset.obs_names.str.split("-").str[0]
subset.obs_names = subset.obs["barcode"].values

# find common barcodes
common = np.intersect1d(subset.obs_names, loom_barcodes)
subset_matched = subset[subset.obs_names.isin(common)].copy() # all matched, 9124

# get cells indices from loom file
# loom_indices = [i for i, bc in enumerate(loom_barcodes) if bc in common]
# sort to ensure matching order
# subset_matched = subset_matched[np.argsort(subset_matched.obs["barcode"])]
# sorted_common = sorted(common)
# loom_indices = [loom_barcodes.index(bc) for bc in sorted_common]
# loom_indices_sorted = sorted(loom_indices)

# make sure loom and h5ad have the same gene symbols
with lp.connect(jj003_loom) as ds:
    loom_genes = ds.ra["Gene"]
# capitalize all gene names in loom to match h5ad
loom_genes = [gene.upper() for gene in loom_genes]

# find common genes
subset_genes = subset_matched.var_names
common_genes = np.intersect1d(loom_genes, subset_genes)

# get indices for common genes and cells
# make sure loom_barcodes and loom_genes are string arrays
loom_barcodes = np.array(loom_barcodes).astype(str)
loom_genes = np.array(loom_genes).astype(str)

# Now build indices safely
loom_cell_idx = []
adata_cell_idx = []

for cell in common:
    matches = np.where(loom_barcodes == cell)[0]
    if matches.size > 0:
        loom_cell_idx.append(matches[0])
        adata_cell_idx.append(subset_matched.obs_names.get_loc(cell))

loom_gene_idx = []
adata_gene_idx = []

for gene in common_genes:
    matches = np.where(loom_genes == gene)[0]
    if matches.size > 0:
        loom_gene_idx.append(matches[0])
        adata_gene_idx.append(subset_matched.var_names.get_loc(gene))

loom_gene_idx_sorted = np.sort(loom_gene_idx)
loom_cell_idx_sorted = np.sort(loom_cell_idx)

# pad missing genes with 0
missing_genes = np.setdiff1d(subset_genes, loom_genes)

with lp.connect(jj003_loom) as ds:
    # slice the data
    spliced_sorted = ds.layers["spliced"][loom_gene_idx_sorted, :][:, loom_cell_idx_sorted]
    unspliced_sorted = ds.layers["unspliced"][loom_gene_idx_sorted, :][:, loom_cell_idx_sorted]
    
    # initialize zero matrices for missing genes
    missing_spliced = np.zeros((len(missing_genes), spliced_sorted.shape[1]))
    missing_unspliced = np.zeros((len(missing_genes), unspliced_sorted.shape[1]))

    # append the missing genes (zero-filled) to the sliced data
    spliced = np.vstack([spliced_sorted, missing_spliced])
    unspliced = np.vstack([unspliced_sorted, missing_unspliced])

    # ensure final shape matches (genes x cells)
    assert spliced.shape[0] == len(subset_genes), "Mismatch in gene count"
    assert unspliced.shape[0] == len(subset_genes), "Mismatch in gene count"

    # assign these matrices to your AnnData object
    subset_matched.layers["spliced"] = spliced.T
    subset_matched.layers["unspliced"] = unspliced.T

# add ambiguous data as well
ambiguous_sorted = np.zeros_like(spliced_sorted)

# pad the ambiguous data with zeros for any missing genes (like for spliced/unspliced)
missing_ambiguous = np.zeros((len(missing_genes), ambiguous_sorted.shape[1]))

# combine the sorted ambiguous data with the missing (zero-filled) data
ambiguous = np.vstack([ambiguous_sorted, missing_ambiguous])

# 4. Ensure the shape matches the expected dimensions (genes x cells)
assert ambiguous.shape[0] == len(subset_genes), "Mismatch in gene count"
assert ambiguous.shape[1] == spliced.shape[1], "Mismatch in cell count"

# 5. Assign the ambiguous layer to your AnnData object
subset_matched.layers["ambiguous"] = ambiguous.T 

#######################################################################################
# INTEGRATING VCY MATRICES INTO MAIN H5AD FILE

adata = sc.read_loom(jj003_loom)
adata1 = sc.read_loom(jj004_loom)
adata2 = sc.read_loom(jj005_loom)
adata_org = sc.read_h5ad(file1)

adata.obs['sample'] = 'JJ003'
adata1.obs['sample'] = 'JJ004'
adata2.obs['sample'] = 'JJ005'

# change cell names in the vly h5ad files
adata.obs.index.name = None
cell_names = adata.obs.index
adata.obs['cellID'] = [c.split(":")[1].rstrip("x") + "-1-0" for c in cell_names]
adata.obs_names = adata.obs['cellID']
adata.obs.index.name = None

adata1.obs.index.name = None
cell_names1 = adata1.obs.index
adata1.obs['cellID'] = [c.split(":")[1].rstrip("x") + "-1-1" for c in cell_names1]
adata1.obs_names = adata1.obs['cellID']
adata1.obs.index.name = None

adata2.obs.index.name = None
cell_names2 = adata2.obs.index
adata2.obs['cellID'] = [c.split(":")[1].rstrip("x") + "-1-2" for c in cell_names2]
adata2.obs_names = adata2.obs['cellID']
adata2.obs.index.name = None

# capitalize all gene names to match that of the main h5ad
adata.var_names = [gene.upper() for gene in adata.var_names]
adata1.var_names = [gene.upper() for gene in adata1.var_names]
adata2.var_names = [gene.upper() for gene in adata2.var_names]

# find shared cells and genes to match that of the main h5ad
# these 10 genes below exist in adata_org but not in adata
#genes_to_change = ['4933427D14RIK-1', 'FAM220A-1', 'GM16364-1', 'GM16499-1',
#                  'GM16701-1', 'GM35438-1', 'GM5089-1', 'GM5468-1', 'PAKAP-1', 'PTP4A1-1']

shared_genes = adata_org.var_names.intersection(adata.var_names) # 18963, should be 18973
shared_cells = adata_org.obs_names.intersection(adata.obs_names) # 9124, all matched
# adata_filtered = adata[shared_cells, shared_genes] -- doesn't work

shared_genes1 = adata_org.var_names.intersection(adata1.var_names) # 18973, all matched
shared_cells_1 = adata_org.obs_names.intersection(adata1.obs_names) # 7250, should be 7264
# adata1_filtered = adata1[shared_cells_1, shared_genes1] -- doesn't work

shared_genes2 = adata_org.var_names.intersection(adata2.var_names) # 18973, all matched
shared_cells_2 = adata_org.obs_names.intersection(adata2.obs_names) # 6870, should be 6911
# adata2_filtered = adata2[shared_cells_2, shared_genes2] -- doesn't work

# integrate into main h5ad, and coercing mismatch cell/gene to 0
def align_layer(source_adata, target_adata, layer_name):
    """
    Aligns a layer from source_adata to the cells and genes in target_adata.
    Missing cells/genes are filled with zeros.
    """
    # get layer from source (fallback to .X if layer is not present)
    X = source_adata.layers[layer_name] if layer_name in source_adata.layers else source_adata.X

    # build a full matrix of zeros
    shape = (target_adata.n_obs, target_adata.n_vars)
    X_aligned = sp.csr_matrix(shape, dtype=X.dtype)

    # find shared genes and cells
    shared_cells = target_adata.obs_names.intersection(source_adata.obs_names)
    shared_genes = target_adata.var_names.intersection(source_adata.var_names)

    # index maps
    target_cell_idx = target_adata.obs_names.get_indexer(shared_cells)
    target_gene_idx = target_adata.var_names.get_indexer(shared_genes)

    source_cell_idx = source_adata.obs_names.get_indexer(shared_cells)
    source_gene_idx = source_adata.var_names.get_indexer(shared_genes)

    # subset layer
    X_sub = X[source_cell_idx[:, None], source_gene_idx]

    # fill aligned matrix
    X_aligned[target_cell_idx[:, None], target_gene_idx] = X_sub

    return X_aligned

# issues with getting index will arise, mainly due to duplicates
for ad in [adata, adata1, adata2, adata_org]:
    ad.var_names_make_unique()

layer_names = ['spliced', 'unspliced', 'ambiguous']

# initialize empty layers
combined_layers = {layer: [] for layer in layer_names}

# loop over each adata and each layer
# there are some duplicates in gene names in adata
adata.var_names_make_unique()

for ad in [adata, adata1, adata2]:
    for layer in layer_names:
        aligned = align_layer(ad, adata_org, layer)
        combined_layers[layer].append(aligned)


# Sum the aligned matrices from all 3 adatas and add to adata_org
for layer in layer_names:
    adata_org.layers[layer] = sum(combined_layers[layer])

adata_org.write_h5ad('/Users/rikac/Documents/project_lab/spring_25/scrna_seq_analysis/bc3306_with_vcy.h5ad')

















































