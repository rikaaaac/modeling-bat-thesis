import scanpy as sc
import anndata as ad
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scrublet as scr
from matplotlib.pyplot import rc_context
import igraph
from scipy.stats import median_abs_deviation
import leidenalg
import openpyxl
import harmonypy
from pandas.plotting import table
sc.logging.print_versions()
sc.settings.verbosity = 3

path = '/Users/rikac/Documents/test/'
matrix_files = '/Users/rikac/Documents/test/matrix_files/'

# COMBINE ALL H5 FILES AND READING THEM

def load_all_samples(h5_path):
    adatas = []

    for filename in os.listdir(h5_path):
        if filename.endswith("_filtered_feature_bc_matrix.h5"):
            samplename = filename.replace("_filtered_feature_bc_matrix.h5", "")
            adata = sc.read_10x_h5(os.path.join(h5_path, filename), gex_only=True)
            adata.var_names_make_unique()
            adata.obs["sample"] = samplename
            adatas.append(adata)

    combined = ad.concat(adatas, join="outer")
    combined.obs_names_make_unique()
    print(combined.obs["sample"].value_counts())
    return combined

adata = load_all_samples(matrix_files)

# QUALITY CONTROL

# mitochondrial genes, 'Mt-' for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

# calculate common QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

# inspect some QC metrics
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    groupby='sample'
)

plt.savefig(path + '/QC/violin_plot_qc.png')
plt.close()

























































