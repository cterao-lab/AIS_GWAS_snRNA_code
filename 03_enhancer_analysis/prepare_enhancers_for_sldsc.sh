#!/bin/bash
set -euo pipefail

INPUT_DIR="results/cluster_active_enhancers"
OUTPUT_DIR="results/cluster_active_enhancers_750bp"

mkdir -p "${OUTPUT_DIR}"

for bed in "${INPUT_DIR}"/*_active_enhancers_CPM1.bed
do
    cluster=$(basename "${bed}" _active_enhancers_CPM1.bed)

    awk 'BEGIN{OFS="\t"}{
        start=$2-750;
        if(start<0) start=0;
        print $1,start,$3+750
    }' "${bed}" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${OUTPUT_DIR}/${cluster}_active_enhancers_CPM1_750bp.bed"
done
