#!/bin/bash
set -euo pipefail

# Build a custom LD reference panel for MAGMA by merging
# 1000 Genomes Phase 3 EAS, EUR, and AMR reference panels.

PLINK=plink

REF_DIR="data/1000G"
OUT_DIR="results/MAGMA_LD_reference"

mkdir -p "${OUT_DIR}"

EAS="${REF_DIR}/g1000_eas"
EUR="${REF_DIR}/g1000_eur"
AMR="${REF_DIR}/g1000_amr"

# ------------------------------------------------------------
# 1. Merge AMR and EUR
# ------------------------------------------------------------

${PLINK} \
  --bfile "${AMR}" \
  --bmerge "${EUR}.bed" "${EUR}.bim" "${EUR}.fam" \
  --make-bed \
  --out "${OUT_DIR}/amr_eur_initial" || true

# Remove variants causing merge conflicts and retry.
if [ -f "${OUT_DIR}/amr_eur_initial-merge.missnp" ]; then

  ${PLINK} \
    --bfile "${AMR}" \
    --exclude "${OUT_DIR}/amr_eur_initial-merge.missnp" \
    --make-bed \
    --out "${OUT_DIR}/amr_filtered"

  ${PLINK} \
    --bfile "${EUR}" \
    --exclude "${OUT_DIR}/amr_eur_initial-merge.missnp" \
    --make-bed \
    --out "${OUT_DIR}/eur_filtered"

  ${PLINK} \
    --bfile "${OUT_DIR}/amr_filtered" \
    --bmerge \
      "${OUT_DIR}/eur_filtered.bed" \
      "${OUT_DIR}/eur_filtered.bim" \
      "${OUT_DIR}/eur_filtered.fam" \
    --make-bed \
    --out "${OUT_DIR}/amr_eur"
else
  mv "${OUT_DIR}/amr_eur_initial.bed" "${OUT_DIR}/amr_eur.bed"
  mv "${OUT_DIR}/amr_eur_initial.bim" "${OUT_DIR}/amr_eur.bim"
  mv "${OUT_DIR}/amr_eur_initial.fam" "${OUT_DIR}/amr_eur.fam"
fi
