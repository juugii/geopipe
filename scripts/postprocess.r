library(Seurat)

project = "GSE135194"

sample <- snakemake@wildcards[["sample"]]
output <- snakemake@output[[1]]

# read gzipped emat
dat <- gzfile('GSE135194.csv.gz')
txt <- readLines(dat)
mat <- read.csv(textConnection(txt), row.names = 1)

seu <- CreateSeuratObject(mat, project = "sample")
# Create a basic processing for UMAP embedding
seu <- NormalizeData(seu)
seu <- ScaleData(seu)
seu <- FindVariableFeatures(seu)
seu <- RunPCA(seu, features = VariableFeatures(seu))
seu <- FindNeighbors(seu)
seu <- FindClusters(seu)
seu <- RunUMAP(seu)
seu <- RunUMAP(seu, dims = 1:10)


# exporting the plot
pdf(output)
	DimPlot(seu, label = TRUE) + NoLegend()
dev.off()

# change with edited snakemake@output[[1]] ?
saveRDS(seu, file = paste0("outs/seu_", sample, ".seurat"))
