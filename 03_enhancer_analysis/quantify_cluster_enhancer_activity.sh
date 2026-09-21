#!/bin/bash
set -euo pipefail

# ReapTEC was run according to Oguchi et al. and the official
# ReapTEC implementation to generate CTSS profiles and a shared
# enhancer catalogue.

ENHANCER_BED="data/reaptec_shared_enhancers.bed"
CLUSTER_CTSS_DIR="data/cluster_CTSS"
OUT_DIR="results/cluster_enhancer_counts"

mkdir -p "${OUT_DIR}"

for ctss_bed in "${CLUSTER_CTSS_DIR}"/*.bed.gz
do
    cluster=$(basename "${ctss_bed}" .bed.gz)

    bedtools coverage \
        -a "${ENHANCER_BED}" \
        -b "${ctss_bed}" \
        -counts \
        > "${OUT_DIR}/${cluster}_enhancer_counts.txt"
done
