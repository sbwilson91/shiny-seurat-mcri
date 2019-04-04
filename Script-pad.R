library(Seurat)
library(tidyverse)
# trial dataset: mafb:gata3 E6 only D7+18 from 2017

seurat <- readRDS("/group/kidn1/Group-Little_MCRI/People/Sean/R/SingleCellAnalysis/Reporters_18vs25Days/output/seurat/D18vsD25_combined_20clusters_seurat.rds")
seurat <- SetIdent(seurat, ident.use = seurat@meta.data$res.1.6)
TSNEPlot(seurat)

markers <- FindAllMarkers(seurat, only.pos = T)
markers %>% filter(cluster == 0)
kidney <- grep(pattern = "[[:punct:]]", x = seurat@data@Dimnames[[1]])

seurat@data@Dimnames[[1]][kidney]


#3d tsne projection
seurat <- RunTSNE(object = seurat, reduction.use = "cca.aligned", dims.use = 1:10, do.fast = TRUE, dim.embed = 3)

tsne_1 <- seurat@dr$tsne@cell.embeddings[,1]

tsne_2 <- seurat@dr$tsne@cell.embeddings[,2]

tsne_3 <- seurat@dr$tsne@cell.embeddings[,3]

library(scatterplot3d)

scatterplot3d(x = tsne_1, y = tsne_2, z = tsne_3, color = as.numeric(1:20)[seurat@ident])

library(rgl) #interactive 3d plotting

plot3d(x = tsne_1, y = tsne_2, z = tsne_3, col = as.numeric(1:10)[seurat@ident], type="s",radius=0.3)

rgl::rglwidget() #save as html