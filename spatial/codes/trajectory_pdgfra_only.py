"""
trajectory analysis using cellrank RealTimeKernel w/ moscot.
runs two analyses:
  1. pdgfra+ cells 
  2. tdtomato+ cells

checkpointing strategy:
  - after tmk.compute_transition_matrix(): write tmat to adata.obsp, save h5ad
  - after g.compute_fate_probabilities(): all gpcca results are in adata, save h5ad
  - to reload: sc.read_h5ad -> RealTimeKernel.from_adata -> GPCCA(tmk)

plotting: wrapped in try/except so failures don't interrupt computation.
"""

import os
# use async cuda allocator to reduce gpu memory fragmentation
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "cuda_async"
# limit jax to 80% of gpu memory so it doesn't try to pre-allocate everything
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.8"
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import geopandas as gpd
from shapely import wkt
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for hpc
import matplotlib.pyplot as plt
import jax
import cellrank as cr
from cellrank.kernels import RealTimeKernel
from moscot.problems.time import TemporalProblem
import harmonypy as hm
import scipy.sparse


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# --- paths -------------------------------------------------------------------

DATA_PATH = "/insomnia001/depts/edu/BIOLBC3141_Fall2025/rc3517/lab/data/combined_samples_analyzed_one_canvas_labeled.h5ad"
OUT_DIR_PDGFRA = "/insomnia001/depts/edu/BIOLBC3141_Fall2025/rc3517/lab/visium/trajectory/pdgfra"
# OUT_DIR_TDTOM = "/insomnia001/depts/edu/BIOLBC3141_Fall2025/rc3517/lab/visium/trajectory/tdtom"
os.makedirs(OUT_DIR_PDGFRA, exist_ok=True)
# os.makedirs(OUT_DIR_TDTOM, exist_ok=True)
PLOT_DIR_PDGFRA = os.path.join(OUT_DIR_PDGFRA, "plots")
# PLOT_DIR_TDTOM = os.path.join(OUT_DIR_TDTOM, "plots")
os.makedirs(PLOT_DIR_PDGFRA, exist_ok=True)
# os.makedirs(PLOT_DIR_TDTOM, exist_ok=True)


# --- helpers -----------------------------------------------------------------

def savefig(path_dir, name):
    path = os.path.join(path_dir, name)
    plt.savefig(path, bbox_inches="tight", dpi=300)
    plt.close()
    log.info("saved plot: %s", path)


def try_plot(fn, plot_dir, name, *args, **kwargs):
    """call fn(*args, **kwargs) and save figure; silently skip on error."""
    try:
        fn(*args, **kwargs)
        savefig(plot_dir, name)
    except Exception as e:
        log.warning("plot '%s' failed, skipping: %s", name, e)
        plt.close()


def checkpoint_tmat(tmk, adata, path):
    """write transition matrix into adata.obsp and save h5ad checkpoint."""
    log.info("writing transition matrix to adata.obsp...")
    tmk.write_to_adata()
    # coarse_fwd stores a nested AnnData which can break older anndata writers;
    # drop it — not needed for fate probabilities or lineage drivers
    adata.uns.pop("coarse_fwd", None)
    if "geometry" in adata.obs.columns:
        adata.obs["geometry"] = adata.obs["geometry"].astype(str)
    adata.write_h5ad(path)
    log.info("checkpoint saved: %s", path)


def checkpoint_gpcca(adata, path):
    """save adata after estimator fit — all results are already in adata."""
    adata.uns.pop("coarse_fwd", None)
    if "geometry" in adata.obs.columns:
        adata.obs["geometry"] = adata.obs["geometry"].astype(str)
    adata.write_h5ad(path)
    log.info("checkpoint saved: %s", path)


def reload_tmk_and_estimator(path):
    """reload tmk and cflare estimator from h5ad checkpoint."""
    log.info("loading checkpoint: %s", path)
    adata = sc.read_h5ad(path)
    tmk = RealTimeKernel.from_adata(adata, key="T_fwd")
    g = cr.estimators.CFLARE(tmk)
    return adata, tmk, g

# --- analysis function -----------------------------------------

def run_trajectory(
    adata,
    cluster_key,
    ba_lineage,
    driver_clusters,
    label,
    eigen_num,
    out_dir,
    plot_dir,
):
    """
    run full trajectory analysis for one adata object.

    args:
        adata:              anndata object (cells x genes), must have 'timepoints' in obs
        cluster_key:        obs column with cell type labels
        ba_lineage:         lineage name to compute drivers for
        driver_clusters:    list of cluster names to restrict driver gene search
        label:              string prefix for log messages and plot filenames
        eigen_num:       # of eigenvalues to compute for cflare
        plot_dir:           directory to save plots
   """

    # marginal scoring 
    log.info("[%s] scoring genes for marginals...", label)
    tp = TemporalProblem(adata)
    tp = tp.score_genes_for_marginals(
        gene_set_proliferation="mouse", gene_set_apoptosis="mouse"
    )

    try_plot(
        lambda: sc.pl.embedding(
            adata, basis="umap",
            color=[cluster_key, "proliferation", "apoptosis"],
            show=False,
        ),
        plot_dir,
        f"{label}_marginals.png",
    )

    # prepare and solve ot problem
    log.info("[%s] preparing temporal problem...", label)
    adata.obs["timepoints_float"] = adata.obs["timepoints"].astype(float)
    # use pca embedding instead of full gene space to reduce memory
    tp = tp.prepare(time_key="timepoints_float", joint_attr="X_pca")

    # solve on cpu to avoid gpu oom on large matrices (27606, 80658)
    log.info("[%s] solving ot problem (rank=500, tau_a=0.95)...", label)
    cpu = jax.devices("cpu")[0]
    with jax.default_device(cpu):
        tp = tp.solve(rank=500, tau_a=0.95, scale_cost="mean")

    # --- build kernel and compute transition matrix --------------------------
    # sparsify on cpu with small batch_size to avoid gpu oom:
    # min_row mode guarantees every row retains at least one entry
    log.info("[%s] building RealTimeKernel (sparsify on cpu)...", label)
    with jax.default_device(cpu):
        tmk = RealTimeKernel.from_moscot(
            tp,
            sparse_mode="min_row",
            sparsify_kwargs={"batch_size": 64},
        )

    log.info("[%s] computing transition matrix...", label)
    tmk.compute_transition_matrix(
        self_transitions="all", conn_weight=0.2, threshold="auto"
    )

    checkpoint_tmat(tmk, adata, os.path.join(out_dir, f"{label}_tmat.h5ad"))

    # visualize kernel
    try_plot(
        lambda: tmk.plot_random_walks(
            max_iter=500,
            start_ixs={"timepoints_float": 11.5},
            basis="X_umap",
            seed=0,
            dpi=150,
            size=30,
            figsize=(10, 7),
            show=False,
        ),
        plot_dir,
        f"{label}_random_walks.png",
    )

    try_plot(
        lambda: tmk.plot_single_flow(
            cluster_key=cluster_key,
            time_key="timepoints_float",
            cluster=ba_lineage,
            min_flow=0.1,
            xticks_step_size=1,
            show=False,
        ),
        plot_dir,
        f"{label}_single_flow.png",
    )

    log.info("[%s] fitting CFLARE (k=10)...", label)
    g = cr.estimators.CFLARE(tmk)
    g.fit(k=eigen_num)

    log.info("[%s] predicting terminal states...", label)
    g.predict(cluster_key=cluster_key, method='kmeans', n_clusters_kmeans=10)

    try_plot(
        lambda: g.plot_macrostates(
            which='terminal', discrete=True, legend_loc="right", s=100, show=False
        ),
        plot_dir,
        f"{label}_terminal_discrete.png",
    )

    try_plot(
        lambda: g.plot_macrostates(
            which='terminal', discrete=False, legend_loc="right", show=False
        ),
        plot_dir,
        f"{label}_terminal_continuous.png",
    )

    # fate probabilities and lineage drivers
    log.info("[%s] computing fate probabilities...", label)
    g.compute_fate_probabilities(tol=1e-9)

    checkpoint_gpcca(adata, os.path.join(out_dir, f"{label}_gpcca.h5ad"))

    all_lineage_names = list(g.fate_probabilities.names)
    matched = [n for n in all_lineage_names if str(n).startswith(ba_lineage)]
    if not matched:
        raise ValueError(
            f"no lineage names starting with {ba_lineage!r}. "
            f"available: {all_lineage_names}"
        )
    log.info("[%s] computing lineage drivers for: %s", label, matched)
    driver_df = g.compute_lineage_drivers(
        lineages=matched,
        cluster_key=cluster_key,
        clusters=driver_clusters,
    )

    out_csv = os.path.join(out_dir, f"{label}_lineage_drivers.csv")
    driver_df.to_csv(out_csv)
    log.info("[%s] lineage drivers saved: %s", label, out_csv)

    return adata, tmk, g, driver_df

# main

def main():
    log.info("loading adata...")
    adata = sc.read_h5ad(DATA_PATH)
    adata.obs["geometry"] = adata.obs["geometry"].apply(wkt.loads)

    mask = adata.obs["pdgfra+"] == True
    adata_pdgfra = adata[mask].copy()
    # remove neurons, heart, neural tube, pharyngeal region
    mask2 = (
        (adata_pdgfra.obs['pdgfra_cell_types'] == 'Heart') |
        (adata_pdgfra.obs['pdgfra_cell_types'] == 'Neural tube') |
        (adata_pdgfra.obs['pdgfra_cell_types'] == 'Pharyngeal region') |
        (adata_pdgfra.obs['pdgfra_cell_types'] == 'Neurons') |
        (adata_pdgfra.obs['pdgfra_cell_types'] == 'Endothelium')
    )
    adata_pdgfra = adata_pdgfra[~mask2].copy()
    log.info("pdgfra subset: %d cells", adata_pdgfra.n_obs)

    # spatial overview plots
    try_plot(
        lambda: gpd.GeoDataFrame(adata.obs, geometry="geometry").plot(
            column="cell_type_labels",
            cmap="tab20",
            ax=plt.subplots(figsize=(10, 8))[1],
            legend=True,
            legend_kwds={"loc": "upper left", "bbox_to_anchor": (1, 1)},
        ),
        PLOT_DIR_PDGFRA,
        "all_spatial_cell_types.png",
    )

    try_plot(
        lambda: sc.pl.umap(adata_pdgfra, show=False),
        PLOT_DIR_PDGFRA,
        "all_umap.png",
    )

    # analysis 1: Pdgfra+ cells
    log.info("=== analysis 1: Pdgfra+ cells ===")

    # also exclude the cell types that were removed earlier
    driver_clusters_all = [
        "Fibroblasts including BA progenitors",
        "Chondrocytes",
        "Endodermal epithelium",
        "Dermal bone",
        "Skeletal muscle",
        "LPM",
        "Epidermis",
        "Dermal fibroblast",
        "Cartilage",
        "Tendon",
    ]

    adata_pdgfra, tmk_all, g_all, _ = run_trajectory(
        adata=adata_pdgfra,
        cluster_key="pdgfra_cell_types",
        ba_lineage="Fibroblasts including BA progenitors",
        driver_clusters=driver_clusters_all,
        label="pdgfra",
        out_dir=OUT_DIR_PDGFRA,
        plot_dir=PLOT_DIR_PDGFRA,
        eigen_num=10,
    )

    log.info("done, output is in %s", OUT_DIR_PDGFRA)

if __name__ == "__main__":
    main()



































