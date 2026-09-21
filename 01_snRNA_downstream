#milo_DA_analysis.R
#-------------------------------
# Differential abundance analysis using miloR (case vs control)
library(miloR)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

# Convert Seurat object to SingleCellExperiment.
# Use RNA assay here to avoid ambiguity; Milo graph construction below uses
# Harmony embeddings rather than integrated expression values.
sce <- as.SingleCellExperiment(seuratObj, assay = "RNA")

reducedDim(sce, "HARMONY") <- seuratObj[["harmony"]]@cell.embeddings
reducedDim(sce, "UMAP") <- seuratObj[["umap"]]@cell.embeddings

# Output directory
output_dir <- "results/milo"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Build Milo object ----
set.seed(13)

milo.obj <- Milo(sce)
milo.obj <- buildGraph(milo.obj, k = 30, d = 30, reduced.dim = "HARMONY")
milo.obj <- makeNhoods(milo.obj, k = 30, d = 30, prop = 0.2, reduced_dims = "HARMONY")

# ---- Count cells per sample ----
meta.data <- seuratObj@meta.data
milo.obj <- countCells(milo.obj, samples = "id", meta.data = meta.data)

# ---- Design matrix ----
design.df <- meta.data %>%
  select(id, status) %>%
  distinct()

design.df <- design.df[match(colnames(milo.obj@nhoodCounts), design.df$id), ]
rownames(design.df) <- design.df$id

design.df$status <- factor(design.df$status, levels = c("ctrl", "case"))

# ---- Differential abundance test ----
milo.obj <- calcNhoodDistance(milo.obj, d = 30, reduced.dim = "HARMONY")

milo.res <- testNhoods(
  milo.obj,
  design = ~ status,
  design.df = design.df
)

# ---- Annotate neighborhoods ----
milo.res <- annotateNhoods(
  milo.obj,
  milo.res,
  coldata_col = "seurat_clusters_annotation"
)

# ---- Save results ----
output_dir <- "results/milo"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  milo.res,
  file = file.path(output_dir, "milo_results.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

# ---- Basic visualization ----
milo.obj <- buildNhoodGraph(milo.obj)

p_umap <- plotNhoodGraphDA(milo.obj, milo.res, layout = "UMAP")

ggsave(
  file.path(output_dir, "milo_umap_DA.pdf"),
  p_umap,
  width = 7,
  height = 6
)

p_beeswarm <- plotDAbeeswarm(
  milo.res,
  group.by = "seurat_clusters_annotation"
)

ggsave(
  file.path(output_dir, "milo_beeswarm_celltype.pdf"),
  p_beeswarm,
  width = 7,
  height = 6
)
#-------------------------------
