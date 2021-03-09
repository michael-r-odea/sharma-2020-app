library(dplyr)
library(Seurat)
library(ggplot2)

sharma <- readRDS("sharma_fourgroups.rds")


#all.markers.wil <- FindAllMarkers(sharma, test.use = "wilcox",
#                                   min.pct = 0, min.cells.feature = 0, min.cells.group = 0, logfc.threshold = 0)

#all.markers <- all.markers.wil


p5m.p0m.all <- FindMarkers(sharma, ident.1 = "P5 - Myelinated", ident.2 = "P0 - Myelinated", test.use = "wilcox", 
                                   min.pct = 0, logfc.threshold = 0, min.cells.group = 0, min.cells.feature = 0)

p5m.p5u.all <- FindMarkers(sharma, ident.1 = "P5 - Myelinated", ident.2 = "P5 - Unmyelinated", test.use = "wilcox", 
                           logfc.threshold = 0, min.pct = 0, min.cells.group = 0, min.cells.feature = 0)

p5m.p0u.all <- FindMarkers(sharma, ident.1 = "P5 - Myelinated", ident.2 = "P0 - Unmyelinated", test.use = "wilcox", 
                           min.pct = 0, logfc.threshold = 0, min.cells.group = 0, min.cells.feature = 0)


p0m.p0u.all <- FindMarkers(sharma, ident.1 = "P0 - Myelinated", ident.2 = "P0 - Unmyelinated", test.use = "wilcox", 
                           min.pct = 0, logfc.threshold = 0, min.cells.group = 0, min.cells.feature = 0)

##

p5u.p0m.all <- FindMarkers(sharma, ident.1 = "P5 - Unmyelinated", ident.2 = "P0 - Myelinated", test.use = "wilcox", 
                           min.pct = 0, logfc.threshold = 0, min.cells.group = 0, min.cells.feature = 0)

p5u.p0u.all <- FindMarkers(sharma, ident.1 = "P5 - Unmyelinated", ident.2 = "P0 - Unmyelinated", test.use = "wilcox", 
                           min.pct = 0, logfc.threshold = 0, min.cells.group = 0, min.cells.feature = 0)


#saveRDS(all.markers, "all_markers.rds") (not using all markers dataset, something weird happens with it)
saveRDS(p5m.p5u.all, "p5m_p5u.rds")
saveRDS(p5m.p0m.all, "p5m_p0m.rds")
saveRDS(p5m.p0u.all, "p5m_p0u.rds")
saveRDS(p0m.p0u.all, "p0m_p0u.rds")
#
saveRDS(p5u.p0m.all, "p5u_p0m.rds")
saveRDS(p5u.p0u.all, "p5u_p0u.rds")

# testing slimmed down seurat object
sharma.slim <- DietSeurat(sharma, counts = FALSE, data = TRUE, scale.data = FALSE, assays = "SCT", dimreducs = "umap")
sharma.slim


# testing that feature/vlnplots are the same with sharma and sharma.slim objects
FeaturePlot(sharma.slim, features = "Cd44")
FeaturePlot(sharma, features = "Cd44")

VlnPlot(sharma.slim, features = "Cd44", pt.size = 0)
VlnPlot(sharma, features = "Cd44", pt.size = 0) + xlab("") + theme(legend.position = "none")

# slimmed object will work
saveRDS(sharma.slim, "sharma_slim.rds")





