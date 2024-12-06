# library
library(Seurat)
library(ggplot2)
library(dplyr)

# setwd
setwd("biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/pub/DRG_spatial/")

# read xenium data
path <- "DRG_Xenium/new_seg"
matching_dirs <- paste0(path, "/drg", c(1:6), "/drg", c(1:6), "/outs/")

xenium_data_list <- list()
for (i in c(1:length(matching_dirs))){
  dir <- paste0("drg", i)
  xenium_path_full <- matching_dirs[i]
  print(xenium_path_full)
  
  xenium.obj <- LoadXenium(xenium_path_full, fov = "fov")
  print(dim(xenium.obj))
  xenium.obj@meta.data$slide_id <- unlist(strsplit(dir, split = "_"))[1]
  
  pdf(paste0("DRG_Xenium/", dir, "_vln_dim.pdf"))
  p <- VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
  print(p)
  p <- ImageDimPlot(xenium.obj, fov = "fov", molecules = c("Gad1", "Sst", "Pvalb", "Gfap"), nmols = 20000)
  print(p)
  p <- ImageFeaturePlot(xenium.obj, fov = "fov", features = "nCount_Xenium", max.cutoff = "q90") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  print(p)
  p <- ImageFeaturePlot(xenium.obj, fov = "fov", features = "nFeature_Xenium", max.cutoff = "q90") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  print(p)
  dev.off()
  
  print(quantile(xenium.obj$nCount_Xenium, seq(0, 1, 0.1)))
  print(quantile(xenium.obj$nFeature_Xenium, seq(0, 1, 0.1)))
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 20)
  print(dim(xenium.obj))
  
  xenium.obj@meta.data$cell_size <- sapply(rownames(xenium.obj@meta.data), FUN = function(x){xenium.obj@images$fov$segmentation[x]@polygons[[1]]@area})
  
  xenium_data_list[[i]] <- xenium.obj
  names(xenium_data_list)[i] <- unlist(strsplit(dir, split = "_"))[1]
  
}

xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- NormalizeData(x)
})

features <- rownames(xenium.obj)
xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "rpca")
xenium.combined <- IntegrateData(anchorset = xenium.anchors)

DefaultAssay(xenium.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
xenium.combined <- ScaleData(xenium.combined, verbose = FALSE)
xenium.combined <- RunPCA(xenium.combined, npcs = 30, verbose = FALSE)
xenium.combined <- RunUMAP(xenium.combined, reduction = "pca", dims = 1:30)
xenium.combined <- FindNeighbors(xenium.combined, reduction = "pca", dims = 1:30)
xenium.combined <- FindClusters(xenium.combined, resolution = 0.3, cluster.name = "xenium_clusters_p3")
xenium.combined <- FindClusters(xenium.combined, resolution = 0.1, cluster.name = "xenium_clusters_p1")

dim(xenium.combined) # > 10: 38462; > 20 :19885; >0: 80269
saveRDS(xenium_data_list, "DRG_Xenium/xenium_data_list_cutoff20.rds")
saveRDS(xenium.combined, "DRG_Xenium/xenium.combined_cutoff20.rds")

xenium.combined <- readRDS("DRG_Xenium/xenium.combined_cutoff20.rds")

pdf("DRG_Xenium/xenium_umap_p1_cutoff20.pdf", width = 6, height = 5)
DimPlot(xenium.combined, label = T)
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3") + ggtitle("Xenium Clusters")
FeaturePlot(xenium.combined, label = T, features = "nCount_Xenium", max.cutoff = "q90")
FeaturePlot(xenium.combined, label = T, features = "nFeature_Xenium")
FeaturePlot(xenium.combined, label = T, features = "cell_size")
dev.off()

# drg_ref markers
drg_ref_markers <- read.csv("DRG_Xenium/comb_type_markers.csv")[, -1]
top_10_vector <- drg_ref_markers %>%
  group_by(cluster) %>%  # Group by columnC
  #arrange(desc(columnB)) %>%  # Arrange in descending order of columnB
  slice_head(n = 5) %>%  # Select the top 10 rows for each group
  pull(gene)

top_10_vector <- unique(top_10_vector)

pdf("DRG_Xenium/top5_dotplot_cutoff20.pdf", width = 11, height = 5)
DotPlot(xenium.combined, features = top_10_vector) + RotatedAxis()
DotPlot(xenium.combined, features = top_10_vector, group.by = "xenium_clusters_p3") + RotatedAxis()
dev.off()

table(xenium.combined$xenium_clusters_p3)
#xenium.combined@meta.data$xenium_clusters_p3_cha <- as.character(xenium.combined@meta.data$xenium_clusters_p3)


xenium.combined@meta.data$xenium_clusters_p3_annot <- recode(xenium.combined@meta.data$xenium_clusters_p3, 
                                                             "0"= "Schwann/SGC", "1"= "Unknown", "2"= "Fibroblast", "3"= "Mrgprd", "4"= "Endothelium/Pericyte",
                                                             "5"= "Neuron_X", "6"= "Neuron_X", "7"= "Neuron_X", "8"= "Th",
                                                             "9"= "Neuron_X", "10"= "Sst")

#xenium.combined@meta.data$xenium_clusters_p3_annot <- recode(xenium.combined@meta.data$xenium_clusters_p3, 
#                                                             "0"= "Schwann/SGC", "1"= "Unknown", "2"= "Mrgprd", "3"= "Neuron_X", "4"= "Fibroblast",
#                                                             "5"= "Endothelium/Pericyte", "6"= "Neuron_X", "7"= "Neuron_X", "8"= "Th",
#                                                             "9"= "Neuron_X", "10"= "Sst", "11"= "Schwann/SGC")

pdf("DRG_Xenium/xenium_umap_p1_cutoff20_annot.pdf", width = 6, height = 4)
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3_annot", repel = T) + ggtitle("Xenium Annotation Using Markers")
dev.off()
pdf("DRG_Xenium/top5_dotplot_cutoff20_annot.pdf", width = 14, height = 4)
DotPlot(xenium.combined, features = top_10_vector, group.by = "xenium_clusters_p3_annot") + RotatedAxis() + ylab("") + xlab("Marker Genes")
dev.off()

saveRDS(xenium.combined, "DRG_Xenium/xenium.combined_cutoff20_annot.rds")
write.csv(xenium.combined@meta.data, "DRG_Xenium/xenium.combined_cutoff20_annot_obs.csv")
saveRDS(xenium.combined@meta.data, "DRG_Xenium/xenium.combined_cutoff20_annot_obs.rds")

#==================================================================================
# handle reference integration
drg_ref <- readRDS("/dfs3b/ruic20_lab/singlecell/tingtiny/DRG_spatial/DRG_comb_ref_neuron1_nonneuron1.rds")
drg_ref <- subset(x = drg_ref, subset = Study != "Deng_2023")
table(drg_ref$Atlas_annotation)
Idents(drg_ref) <- "Atlas_annotation"
drg_ref <- subset(x = drg_ref, downsample = 2000)
table(drg_ref$Atlas_annotation)
dim(drg_ref)

Idents(drg_ref) <- "Study"
# split the dataset into a list of two seurat objects (stim and CTRL)
drg_ref.list <- SplitObject(drg_ref, split.by = "Study")

# normalize and identify variable features for each dataset independently
drg_ref.list <- lapply(X = drg_ref.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- rownames(xenium.combined)
features <- intersect(features, rownames(drg_ref))
features <- c(features, SelectIntegrationFeatures(object.list = drg_ref.list))
features <- unique(features)

drg_ref.list <- lapply(X = drg_ref.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

drg_ref.anchors <- FindIntegrationAnchors(object.list = drg_ref.list, anchor.features = features, reduction = "rpca")
drg_ref.combined <- IntegrateData(anchorset = drg_ref.anchors)

DefaultAssay(drg_ref.combined) <- "integrated"
drg_ref.combined <- ScaleData(drg_ref.combined, verbose = FALSE)
drg_ref.combined <- RunPCA(drg_ref.combined, npcs = 30, verbose = FALSE)
drg_ref.combined <- RunUMAP(drg_ref.combined, reduction = "pca", dims = 1:30)
drg_ref.combined <- FindNeighbors(drg_ref.combined, reduction = "pca", dims = 1:30)
drg_ref.combined <- FindClusters(drg_ref.combined, resolution = 0.5)

saveRDS(drg_ref.combined, "DRG_Xenium/drg_ref.combined_2k_plus_panel.rds")

pdf("DRG_Xenium/drg_ref_umaps_rpca_2k_plus_panel.pdf", width = 10, height = 5)
DimPlot(drg_ref.combined)
DimPlot(drg_ref.combined, group.by = "Atlas_annotation")
DimPlot(drg_ref.combined, group.by = "Atlas_annotation", label = T)
DimPlot(drg_ref.combined, group.by = "Study", label = T)
dev.off()

#==================================================================================
# co-embedding
#xenium.combined <- readRDS("DRG_Xenium/xenium.combined.rds")
drg_ref.combined <- readRDS("DRG_Xenium/drg_ref.combined_2k_plus_panel.rds")

drg_ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(drg_ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
DefaultAssay(drg_ref.combined)
DefaultAssay(xenium.combined)

obj_list <- list("ref" = drg_ref.combined, "xenium" = xenium.combined)
#obj_list <- lapply(X = obj_list, FUN = function(x) {
#  x <- NormalizeData(x)
#})
features <- rownames(xenium.combined)
features <- intersect(features, rownames(drg_ref.combined))
#obj_list <- lapply(X = obj_list, FUN = function(x) {
#  x <- ScaleData(x, features = features, verbose = FALSE)
#  x <- RunPCA(x, features = features, verbose = FALSE)
#})

coembed.anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = "rpca")
coembed.combined <- IntegrateData(anchorset = coembed.anchors, new.assay.name = "coembed")
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(coembed.combined) <- "coembed"

# Run the standard workflow for visualization and clustering
coembed.combined <- ScaleData(coembed.combined, verbose = FALSE)
coembed.combined <- RunPCA(coembed.combined, npcs = 30, verbose = FALSE)
coembed.combined <- RunUMAP(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindNeighbors(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindClusters(coembed.combined, resolution = 0.5, cluster.name = "Rp5_clusters")
#coembed.combined <- FindClusters(coembed.combined, resolution = 0.8, cluster.name = "Rp8_clusters")

saveRDS(coembed.combined, "DRG_Xenium/coembed.combined_inte.rds")

