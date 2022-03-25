# get metadata from GEO datasets
# Jules GILET <jules.gilet@inserm;fr>

library(GEOquery)

project = "GSE135194"

geo <- getGEO(project)
tb <- pData(geo[[1]])
grep("title|subject|Sex", colnames(tb))
meta <- tb[ , grep("title|subject|Sex|relation", colnames(tb))]

write.csv(meta, file = meta.csv, row.names = FALSE)
