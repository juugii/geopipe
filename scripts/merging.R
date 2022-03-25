# Default Datasets Integration with Seurat v4
# Jules GILET jules.gilet@inserm.fr<>

input <- snakemake@input[["list"]]
gse <- snakemake@params[['gse']]
outz <- snakemake@output[["file"]]
figz <- snakemake@output[["plot"]]
logz <- snakemake@output[["log"]]
countz <- snakemake@output[["matrix"]]

# outz = paste0("outs/seurat/integrated/seu_int_", gse, ".seurat")
# figz = paste0("outs/seurat/integrated/UMAPs_", gse, ".pdf")
# logz = paste0("outs/seurat/integrated/integration_", gse, ".log")

fcon <- file(logz)
sink(fcon, type = "output", append = TRUE)
sink(fcon, type = "message", append = TRUE)

library(Seurat)
library(dplyr)

# reading the induividual samples
# already preprocessed
seu_list <- list()
i <- 1
seu_list[[i]] <- readRDS(input[1])
for (i in 2:length(input)){
	seu_list[[i]] <- readRDS(input[i])
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
pdf(figz)
	DimPlot(seu_int, label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "Run", label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "Age", label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "sex", label = TRUE) + NoLegend()
	DimPlot(seu_int, group.by = "subject", label = TRUE) + NoLegend()
dev.off()

saveRDS(seu_int, outz)

# export read count matrix
mat <- GetAssayData(seu_int, slot = 'counts', assay = "RNA")
write.csv(as.matrix(mat), file = countz)

print("Done.")
sessionInfo()
