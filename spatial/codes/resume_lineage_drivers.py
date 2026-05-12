"""
resume lineage driver computation from pdgfra_gpcca.h5ad checkpoint.
run this instead of rerunning the full script.
"""

import os
import logging
import scanpy as sc
import cellrank as cr
from cellrank.kernels import RealTimeKernel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

CHECKPOINT = "./pdgfra_new/pdgfra_gpcca.h5ad"
OUT_DIR = "./pdgfra_new"
CLUSTER_KEY = "pdgfra_cell_types"
BA_LINEAGE = "Fibroblasts including BA progenitors"

# cell types that were excluded from adata_pdgfra via mask2
DRIVER_CLUSTERS = [
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

log.info("loading checkpoint: %s", CHECKPOINT)
adata = sc.read_h5ad(CHECKPOINT)

tmk = RealTimeKernel.from_adata(adata, key="T_fwd")
g = cr.estimators.CFLARE(tmk)

# fate probabilities are already stored in adata.obsm['lineages_fwd']
# — just need to attach them back to the estimator
g._read_fate_probabilities(adata)

all_lineage_names = list(g.fate_probabilities.names)
log.info("available lineages: %s", all_lineage_names)

matched = [n for n in all_lineage_names if str(n).startswith(BA_LINEAGE)]
if not matched:
    raise ValueError(
        f"no lineage names starting with {BA_LINEAGE!r}. "
        f"available: {all_lineage_names}"
    )

log.info("computing lineage drivers for: %s", matched)
driver_df = g.compute_lineage_drivers(
    lineages=matched,
    cluster_key=CLUSTER_KEY,
    clusters=DRIVER_CLUSTERS,
)

out_csv = os.path.join(OUT_DIR, "pdgfra_lineage_drivers.csv")
driver_df.to_csv(out_csv)
log.info("lineage drivers saved: %s", out_csv)
