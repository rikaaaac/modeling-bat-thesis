import scanpy as sc
import cellrank as cr
from cellrank.kernels import RealTimeKernel
from moscot.problems.time import TemporalProblem
import matplotlib.pyplot as plt
from shapely import wkt
import geopandas as gpd
import scanpy.external as sce


# --- load data ---

adata = sc.read('../../combined_samples_analyzed_one_canvas_labeled.h5ad')

adata.obs['geometry'] = adata.obs['geometry'].apply(wkt.loads)
shapes = gpd.GeoDataFrame(adata.obs, geometry='geometry')


# --- overview plots ---

fig, ax = plt.subplots(figsize=(10, 8))
shapes.plot(column='timepoints', cmap='tab20', ax=ax, legend=True,
            legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

fig, ax = plt.subplots(figsize=(14, 12))
shapes.plot(column='cell_type_labels', cmap='tab20', ax=ax, legend=True,
            legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

sc.pl.umap(adata, color=['timepoints', 'Pdgfra', 'cell_type_labels'],
           frameon=False, ncols=2, color_map='Blues', vmax=3)


# --- marker gene expression ---

marker_genes = [
    'Krt5', 'Krt14', 'Cdh5', 'Kdr', 'Osr2', 'Crabp1', 'Twist2', 'Acta2',
    'Tagln', 'Col2a1', 'Cldn11', 'Foxc1', 'Scx', 'Tnmd', 'Pparg', 'Cebpa',
    'Pecam1', 'Pdgfra',
]
sc.pl.dotplot(adata, marker_genes, groupby='cell_type_labels', standard_scale='var')

# epidermis markers: Krt5, Krt14
for g in ['Krt5', 'Krt14']:
    gene_idx = adata.var_names.get_loc(g)
    shapes['expr'] = adata.X[:, gene_idx].toarray().flatten()
    fig, ax = plt.subplots(figsize=(10, 8))
    shapes.plot(column='expr', cmap='Reds', ax=ax, legend=True, vmax=1)
    ax.set_title(g)
    plt.show()

# endothelium markers: Cdh5, Kdr
for g in ['Kdr', 'Cdh5']:
    gene_idx = adata.var_names.get_loc(g)
    shapes['expr'] = adata.X[:, gene_idx].toarray().flatten()
    fig, ax = plt.subplots(figsize=(10, 8))
    shapes.plot(column='expr', cmap='Reds', ax=ax, legend=True, vmax=1)
    ax.set_title(g)
    plt.show()

# MCT markers
for g in ['Osr2', 'Crabp1']:
    gene_idx = adata.var_names.get_loc(g)
    shapes['expr'] = adata.X[:, gene_idx].toarray().flatten()
    fig, ax = plt.subplots(figsize=(10, 8))
    shapes.plot(column='expr', cmap='Reds', ax=ax, legend=True, vmax=1)
    ax.set_title(g)
    plt.show()

# dermis markers: Twist2
g = 'Twist2'
gene_idx = adata.var_names.get_loc(g)
shapes['expr'] = adata.X[:, gene_idx].toarray().flatten()
fig, ax = plt.subplots(figsize=(10, 8))
shapes.plot(column='expr', cmap='Reds', ax=ax, legend=True, vmax=1)
ax.set_title(g)
plt.show()


# --- pdgfra+ cells ---

mask_pdgfra = adata.obs['pdgfra+'] == True
pdgfra = adata[mask_pdgfra].copy()

shapes_pdgfra = gpd.GeoDataFrame(pdgfra.obs, geometry='geometry')

fig, ax = plt.subplots(figsize=(110, 8))
shapes_pdgfra.plot(column='cell_type_labels', cmap='tab20', ax=ax, legend=True,
                   legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

# original umap before recalculating coordinates
sc.pl.umap(pdgfra, color=['timepoints', 'Pdgfra', 'pdgfra_cell_types'], ncols=2, vmax=3)

# rerun pca
sc.pp.pca(pdgfra, n_comps=50)
# rerun harmony
sce.pp.harmony_integrate(pdgfra, key="sample", basis="X_pca", max_iter_harmony=20)
# rerun neighbors + umap
sc.pp.neighbors(pdgfra, use_rep="X_pca_harmony")
sc.tl.umap(pdgfra, random_state=0)

sc.pl.umap(pdgfra, color=['timepoints', 'pdgfra_cell_types', 'sample'], ncols=2, frameon=False)

sc.tl.rank_genes_groups(pdgfra, groupby='pdgfra_cell_types', use_raw=False,
                        method='wilcoxon', key_added='pdgfra_degs')
sc.pl.rank_genes_groups(pdgfra, key='pdgfra_degs')

sc.pl.umap(pdgfra, color=['timepoints', 'pdgfra_cell_types'])
sc.pl.umap(pdgfra, color=['Pparg', 'Cebpa', 'Gata6', 'Cdh4', 'Sox9', 'Col2a1'],
           vmax=1, color_map='Blues')
           


# --- pdgfra+ (pick out clusters) ---

mask2 = (
    (pdgfra.obs['cell_type_labels'] == 'Cartilage') |
    (pdgfra.obs['cell_type_labels'] == 'LPM') |
    (pdgfra.obs['cell_type_labels'] == 'Tendon') |
    (pdgfra.obs['cell_type_labels'] == 'Dermal fibroblast') |
    (pdgfra.obs['cell_type_labels'] == 'Fibroblasts including BA progenitors') |
    (pdgfra.obs['cell_type_labels'] == 'Chondrocytes') |
    (pdgfra.obs['cell_type_labels'] == 'Dermal bone')
)
adata_pdgfra_clusters = pdgfra[mask2].copy()

shapes_pdgfra_clusters = gpd.GeoDataFrame(adata_pdgfra_clusters.obs, geometry='geometry')

sc.pl.umap(adata_pdgfra_clusters, color=['Pdgfra', 'timepoints', 'cell_type_labels'],
           ncols=2, vmax=3, color_map='Blues', frameon=False)

fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
shapes_pdgfra_clusters.plot(column='cell_type_labels', cmap='tab20', ax=ax, legend=True,
                   legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

# recalculate coordinates
sc.tl.pca(adata_pdgfra_clusters, n_comps=30)
# rerun harmony
sce.pp.harmony_integrate(adata_pdgfra_clusters, key="sample", basis="X_pca", max_iter_harmony=20)
sc.pp.neighbors(adata_pdgfra_clusters, use_rep='X_pca_harmony')
sc.tl.umap(adata_pdgfra_clusters, random_state=42)

sc.pl.umap(adata_pdgfra_clusters, color=['Pdgfra', 'timepoints', 'cell_type_labels'],
           color_map='Blues', vmax=3, ncols=2, frameon=False)

# # find degs
# sc.tl.rank_genes_groups(adata_pdgfra_clusters, groupby='cell_type_labels', 
#                         use_raw=False, layer='norm_log', 
#                         key_added='pdgfra_clusters_deg', method='wilcoxon')
# sc.pl.rank_genes_groups(adata_pdgfra_clusters, key='pdgfra_clusters_deg', n_genes=25)
# # save to csv
# result = sc.get.rank_genes_groups_df(adata_pdgfra_clusters, group=None, key='pdgfra_clusters_deg')
# result.to_csv('pdgfra_clusters_deg.csv', index=False)

# reclustering (0.6)
sc.tl.leiden(adata_pdgfra_clusters, resolution=0.6, key_added='pdgfra_leiden_0.6')

sc.pl.umap(adata_pdgfra_clusters, color='pdgfra_leiden_0.6', frameon=False)

fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
shapes_pdgfra_clusters.plot(column='pdgfra_leiden_0.6', cmap='tab20', ax=ax, legend=True,
                   legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

# find degs for leiden 0.6
sc.tl.rank_genes_groups(adata_pdgfra_clusters, groupby='pdgfra_leiden_0.6', 
                        use_raw=False, layer='norm_log', 
                        key_added='pdgfra_leiden_0.6_deg', method='wilcoxon')
result2 = sc.get.rank_genes_groups_df(adata_pdgfra_clusters, group=None, key='pdgfra_leiden_0.6_deg')
result2.to_csv('../../../pdgfra_leiden_0.6_deg.csv', index=False)

# reclustering (0.7, 0.8)
sc.tl.leiden(adata_pdgfra_clusters, resolution=0.7, key_added='pdgfra_leiden_0.7')
sc.tl.leiden(adata_pdgfra_clusters, resolution=0.8, key_added='pdgfra_leiden_0.8')
sc.pl.umap(adata_pdgfra_clusters, color=['pdgfra_leiden_0.7', 
                                        'pdgfra_leiden_0.8'], frameon=False)

# cluster identities for leiden_0.6
clusters_leiden_06 = {
    '0': 'Fibroblasts and BA progenitors',
    '1': 'Chondrocytes',
    '2': 'Dermal fibroblasts',
    '3': 'E12.5 LPM 1',
    '4': 'Tendon',
    '5': 'E12.5 LPM 2',
    '6': 'E12.5 LPM 3',
    '7': 'Fibroblasts',
    '8': 'E11.5 LPM 1',
    '9': 'E11.5 LPM 2',
    '10': 'Dermal bone'
}

adata_pdgfra_clusters.obs['pdgfra_leiden_0.6_cell_types'] = adata_pdgfra_clusters.obs['pdgfra_leiden_0.6'].map(clusters_leiden_06)

fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
shapes_pdgfra_clusters.plot(column='pdgfra_leiden_0.7', cmap='tab20', ax=ax, legend=True,
                   legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
shapes_pdgfra_clusters.plot(column='pdgfra_leiden_0.8', cmap='tab20', ax=ax, legend=True,
                   legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

# find degs for leiden 0.7, 0.8
sc.tl.rank_genes_groups(adata_pdgfra_clusters, groupby='pdgfra_leiden_0.7', 
                        use_raw=False, layer='norm_log', 
                        key_added='pdgfra_leiden_0.7_deg', method='wilcoxon')
result2 = sc.get.rank_genes_groups_df(adata_pdgfra_clusters, group=None, key='pdgfra_leiden_0.7_deg')
result2.to_csv('../../../pdgfra_leiden_0.7_deg.csv', index=False)

sc.tl.rank_genes_groups(adata_pdgfra_clusters, groupby='pdgfra_leiden_0.8', 
                        use_raw=False, layer='norm_log', 
                        key_added='pdgfra_leiden_0.8_deg', method='wilcoxon')
result2 = sc.get.rank_genes_groups_df(adata_pdgfra_clusters, group=None, key='pdgfra_leiden_0.8_deg')
result2.to_csv('../../../pdgfra_leiden_0.8_deg.csv', index=False)












# --- tdtom cells ---

mask_tdtom = adata.obs['tdtomato+'] == True
tdtom = adata[mask_tdtom].copy()

shapes_tdtom = gpd.GeoDataFrame(tdtom.obs, geometry='geometry')

# new cluster mapping for tdtom leiden 0.6
clusters_mapping = {
    '0': 'Mesenchyme including BA',
    '1': 'E13.5 Dermis',
    '2': 'E12.5-E13.5 Mesenchyme',
    '3': 'E13.5 Pharyngeal region',
    '4': 'E13.5 Pharyngeal region',
    '5': 'E13.5 Skeletal muscle',
    '6': 'E13.5 Dermal bone',
    '7': 'E13.5 Cartilage and tendon',
    '8': 'E12.5 Skeletal muscle',
    '9': 'E11.5 Myotome and LPM',
    '10': 'E13.5 Endoderm',
    '11': 'E13.5 Ganglia',
    '12': 'Meninges',
    '13': 'E12.5 Pharyngeal region',
    '14': 'E12.5 Dorsal Dermis',
    '15': 'Neurons',
    '16': 'Heart',
    '17': 'Neural tube',
    '18': 'Neural tube',
    '19': 'Neural tube',
    '20': 'Smooth muscle',
    '21': 'Other mesenchyme'
}

tdtom.obs['tdtom_leiden_0.6_cell_types'] = tdtom.obs['tdtom_reclustered_leiden_0.6'].map(clusters_mapping)


# remove some cell types from tdtom
mask3 = (
    (tdtom.obs['tdtom_leiden_0.6_cell_types'] == 'Heart') |
    (tdtom.obs['tdtom_leiden_0.6_cell_types'] == 'E13.5 Pharyngeal region') |
    (tdtom.obs['tdtom_leiden_0.6_cell_types'] == 'Neural tube') |
    (tdtom.obs['tdtom_leiden_0.6_cell_types'] == 'Neurons') |
    (tdtom.obs['tdtom_leiden_0.6_cell_types'] == 'E12.5 Pharyngeal region') |
    (tdtom.obs['tdtom_leiden_0.6_cell_types'] == 'E13.5 Ganglia') |
    (tdtom.obs['tdtom_leiden_0.6_cell_types'] == 'E13.5 Endoderm')
)
adata_tdtom = tdtom[~mask3].copy()

# recalculate coordinates
# rerun pca
sc.pp.pca(adata_tdtom, n_comps=30)
# rerun harmony
sce.pp.harmony_integrate(adata_tdtom, key="sample", basis="X_pca", max_iter_harmony=20)
# rerun neighbors + umap
sc.pp.neighbors(adata_tdtom, use_rep="X_pca_harmony")
sc.tl.umap(adata_tdtom, random_state=0)

sc.pl.umap(adata_tdtom, color=['timepoints', 'tdtom_leiden_0.6_cell_types', 'sample'], ncols=2, frameon=False)

# convert geometry to WKT strings before saving
if 'geometry' in adata_tdtom.obs.columns:
    adata_tdtom.obs['geometry'] = adata_tdtom.obs['geometry'].astype(str)
adata_tdtom.write_h5ad('../../tdtomato_labeled.h5ad', compression='gzip')














