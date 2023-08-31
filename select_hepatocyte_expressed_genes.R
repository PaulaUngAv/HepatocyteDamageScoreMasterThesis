library(Seurat)
library(sctransform)
library(ggplot2)

# need to clarify where input file comes from 
lca.mouseMerged <- readRDS(file = 'scRNAseq/Liver Cell Atlas Data/merged/countTable_mouseMergedConditions_randomsubset.rds')

annotation.mouseStSt.ForMerging <- 
  read.csv(file = 'scRNAseq/Liver Cell Atlas Data/merged/annot_mouseStSt_ForMerging_randomsubset.csv')

annotation.mouseNafld.ForMerging <- 
  read.csv(file = 'scRNAseq/Liver Cell Atlas Data/merged/annot_mouseNafld_ForMerging_randomsubset.csv')

## add cell type identities from annotations to seurat object
# if true:
check_identities <- gsub("StSt_", "", names(lca.mouseMerged$orig.ident))
check_identities <- gsub("Nafld_", "", check_identities)

identical(x = annotation.mouseStSt.ForMerging$cell , 
          y = check_identities[1:25000])

identical(x = annotation.mouseNafld.ForMerging$cell , 
          y = check_identities[25001:50000])
new_idents 
# then:

lca.mouseMerged@meta.data$cell_type <- c(annotation.mouseStSt.ForMerging$annot,
                                         annotation.mouseNafld.ForMerging$annot)



#Quality control 
lca.mouseMerged[["percent.mt"]] <- PercentageFeatureSet(lca.mouseMerged, 
                                                             pattern = "^mt-")
VlnPlot(lca.mouseMerged, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt" ), 
        ncol = 3,pt.size = F) 

VlnPlot(lca.mouseMerged, group.by = 'cell_type', features = "nFeature_RNA", pt.size = F ) + 
  theme(legend.position = 'none')
#upper limit 3000-3500, lower limit 200

VlnPlot(lca.mouseMerged, features = "nCount_RNA", group.by = 'cell_type', pt.size = F) + 
  theme(legend.position = 'none')
# 15 000 upper limit

VlnPlot(lca.mouseMerged, features = "percent.mt", group.by = 'cell_type', pt.size = F ) + 
  theme(legend.position = 'none')
#upper limit 7.5

## Choose filtering parameters
# nFeature_RNA > ... & nFeature_RNA < ... & percent.mt < ... & nCount_RNA < ...
## filter in the cluster and sct transform and find neighbors and clusters
# also dimension reduction
# reload data 

lca.mouseMerged <-  readRDS(file = 'scRNAseq/Liver Cell Atlas Data/countTable_mouseMergedConditions_randomsubset_sctransformed.rds')

# plot dimension redution control 
# title?? 
DimPlot(lca.mouseMerged, reduction = 'umap',  
        label = T, group.by = 'seurat_clusters' )

DimPlot(lca.mouseMerged, reduction = 'pca')


## calculate parameters for gene filtering

# 1. Detection Rate in hepatocytes and in other cells and Log Fold Change 
other_cell_types <- unique(annotation.mouseStSt.ForMerging$annot)
unique(annotation.mouseNafld.ForMerging$annot)
# add the ones missing 
other_cell_types <-append()
#remove hepatocytes from list
other_cell_types <- other_cell_types[-5]

lca.mouseMerged.LFC.DetectionRate <- FoldChange(
  lca.mouseMerged,
  ident.1 = 'Hepatocytes',
  ident.2 = other_cell_types,
  assay = 'SCT',
  base = 2,
)
# 2. Average Expression  per cell type and avg. across all cell types but hep.
lca.mouseMerged.averageExpression <- AverageExpression(object = lca.mouseMerged,
                                                            assays = 'SCT')

lca.mouseMerged.averageExpression <- as.data.frame(lca.mouseMerged.averageExpression$SCT)

summary(lca.mouseMerged.averageExpression)

plot(lca.mouseMerged.averageExpression[,1],lca.mouseMerged.Subset.averageExpression[,2])

# average expression across all cell types except hepatocytes (col 5)???
lca.mouseMerged.averageExpression$global_avg_minus_hep <- 
  rowSums(lca.mouseMerged.averageExpression[,-5])

## Universal Hepatocyte Damage Signature Unfiltered 
HDS_unfiltered <- read.csv(file = 'universal_damage_signature_16.03.22.csv',
                              sep = ';')

# Filtering Illustration Plots 


# Filtering 
HDS_filtered 

