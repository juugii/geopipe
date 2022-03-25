library(Seurat)
library(dplyr)

files <- snakemake@input[[1]]
output <- snakemake@output[[1]]
gse <- snakemake@params[['gse']]

# reading the induividual samples
# already preprocessed
seu_list <- list()
seu_list <- readRDS(files[[1]])
for (i in 2:length(files)){
	seu_list[[i]] <- readRDS(files[[i]])
}

# select features
features <- SelectIntegrationFeatures(seu_list)
# identify anchors
anchors <- FindIntegrationAnchors(seu_list, anchor.features = features)

# build integrated dataset
seu_int <- IntegrateData(anchors)
DefaultAssay(seu_int) <- "integrated"

# Postprocessing
seu_int <- ScaleData(seu_int)
seu_int <- RunPCA(seu_int)
seu_int <- RunUMAP(seu_int, reduction = "pca", dims = 1:10)
seu_int <- FindNeighbors(seu_int, reduction = "pca", dims = 1:10)
seu_int <- FindClusters(seu_int)

# Ouputs
pdf(paste0("UMAP_", gse, ".pdf"))
	DimPlot(seu_int, label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "Run", label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "Age", label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "sex", label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "subject", label = TRUE) + NoLegend()
dev.off()

saveRDS(seu_int, file = paste0("outs/seu_int_", gse, ".rds"))
