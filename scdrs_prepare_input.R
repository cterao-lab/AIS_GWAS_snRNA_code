library(Seurat)
library(SeuratDisk)

# -----------------------------
# Prepare Seurat object for scDRS
# -----------------------------

DefaultAssay(seuratObj) <- "RNA"

# Merge Seurat v5 count layers across samples
seuratObj <- JoinLayers(
  seuratObj,
  assay = "RNA"
)

# Extract raw counts
counts_mat <- GetAssayData(
  seuratObj,
  assay = "RNA",
  layer = "counts"
)

# Rebuild a clean RNA assay containing raw counts
seuratObj[["RNA_clean"]] <- CreateAssayObject(
  counts = counts_mat
)

DefaultAssay(seuratObj) <- "RNA_clean"

seuratObj[["RNA"]] <- NULL

seuratObj <- RenameAssays(
  seuratObj,
  RNA_clean = "RNA"
)

# -----------------------------
# Export to h5ad
# -----------------------------

SaveH5Seurat(
  seuratObj,
  filename = "data/seuratObj.h5Seurat",
  overwrite = TRUE
)

Convert(
  "data/seuratObj.h5Seurat",
  dest = "h5ad",
  overwrite = TRUE
)

# -----------------------------
# Export covariates for scDRS
# -----------------------------

meta <- seuratObj@meta.data

covariates <- data.frame(
  nuclei_id = rownames(meta),
  sample = meta$id,
  nFeature = meta$nFeature_RNA
)

write.table(
  covariates,
  file = "data/seuratObj.cov",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
