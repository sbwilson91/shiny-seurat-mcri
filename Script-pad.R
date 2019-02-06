library(Seurat)
library(tidyverse)
# trial dataset: mafb:gata3 E6 only D7+18 from 2017

seurat <- readRDS("E:/R/180428_SingleCellCombination/output/seurat/E6_Seurat.Rds")

TSNEPlot(seurat)

markers <- FindAllMarkers(seurat, only.pos = T)
markers %>% filter(cluster == 0)
kidney <- grep(pattern = "[[:punct:]]", x = seurat@data@Dimnames[[1]])

seurat@data@Dimnames[[1]][kidney]
