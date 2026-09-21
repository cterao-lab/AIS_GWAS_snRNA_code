# Subclustering of FAPs 
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

set.seed(13)

# -------------------------------
# Subset FAPs 
# -------------------------------

faps <- subset(
  seuratObj,
  subset = seurat_clusters_annotation %in% c("FAPs")
)

DefaultAssay(faps) <- "RNA"

# -------------------------------
# Reprocessing and subclustering
# -------------------------------

faps <- NormalizeData(faps)
faps <- FindVariableFeatures(faps, nfeatures = 2000)
faps <- ScaleData(faps)
faps <- RunPCA(
  faps,
  features = VariableFeatures(faps)
)

# Correct donor-level effects
faps <- RunHarmony(
  faps,
  group.by.vars = "id"
)

faps <- RunUMAP(faps, reduction = "harmony", dims=1:20)

#clustering
faps <- FindNeighbors(faps, reduction = "harmony", dims=1:20)
faps <- FindClusters(faps, resolution = 0.1) 

# -------------------------------
# Save output
# -------------------------------

output_dir <- "results/FAP_subclustering"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  faps,
  file = file.path(output_dir, "FAP_subclustered_seurat_object.rds")
)
