# Individual Simple Preprocessing with Seurat
# Jules GILET jules.gilet@inserm.fr<>

input <- snakemake@input[[1]]
output <- snakemake@output[["file"]]
fig <- snakemake@output[["plot"]]
log <- snakemake@output[["log"]]
sample <- snakemake@wildcards[["sample"]]

fcon <- file(log)
sink(fcon, type = "output", append = TRUE)
sink(fcon, type = "message", append = TRUE)


library(Seurat)
library(GEOquery)
library(SRAdb)

# read h5
# for some SRR h5 is missing...
mat <- Read10X_h5(input)

# reading metadata
meta <- read.csv('meta.csv')[ , c("Run", "Age", "sex", "subject")]
meta <- meta[ meta$Run == sample, ]
meta <- meta[ rep(1, ncol(mat)), ]
rownames(meta) <- colnames(mat)

seu <- CreateSeuratObject(mat, project = sample, meta.data = meta)

# Create a basic processing for UMAP embedding
seu <- NormalizeData(seu)
seu <- ScaleData(seu)
seu <- FindVariableFeatures(seu)
seu <- RunPCA(seu, features = VariableFeatures(seu))
seu <- FindNeighbors(seu)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:10)


# exporting the plot
pdf(fig)
	DimPlot(seu, label = TRUE) + NoLegend()
dev.off()

# change with edited snakemake@output[[1]] ?
saveRDS(seu, file = output)

print("Done.")
sessionInfo()
