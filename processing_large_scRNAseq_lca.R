### Filter Universal Damage Signature Gene List numerically
# script to run on cluster 
library(Seurat)
library(sctransform)
library(ggplot2)

path_data =  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/'
path_annotation = path_data =  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseStStAll.csv'

#annotation.mouseStSt <- read.csv(
#  file = 'scRNAseq/Liver Cell Atlas Data /annot_mouseStStAll.csv')

annotation.mouseStSt <- read.csv(
  file = path_annotation)

counts <- ReadMtx(mtx =  paste0(path_data,'matrix.mtx.gz', sep = '' ), 
                  cells = paste0(path_data,'barcodes.tsv.gz', sep = '' ),
                  features = paste0(path_data,'features.tsv.gz', sep = '' ),
                  feature.column = 1) 


# plot 

## filtering 

## a) create subset small enough to 'play with' 
# using equivalent of sample function of Seutat -> randomly select 50 000 cells 

sample_indices <-sample(c(1:length(annotation.mouseStSt$cell)), 50000, 
                        replace = F )

annotation.mouseStSt.Subset <-  annotation.mouseStSt[sample_indices, ]

lca.mouseStSt <- CreateSeuratObject(counts = counts, 
                                    project = "liver_cell_atlas_stst", 
                                    min.cells = 3)

lca.mouseStSt.Subset <- subset(lca.mouseStSt, 
                               cells = annotation.mouseStSt.Subset[,'cell'] )

write.csv(annotation.mouseStSt.Subset, file = 'annot_mouseStStAll_randomsubset.csv')
saveRDS(lca.mouseStSt.Subset, file = 'countTable_mouseStSt_randomsubset.rds')

##########################################################
# work with subset 

lca.mouseStSt.Subset <- readRDS(file = 'scRNAseq/Liver Cell Atlas Data/countTable_mouseStSt_randomsubset.rds')
lca.mouseNafld.Subset <- readRDS(file = 'scRNAseq/Liver Cell Atlas Data/countTable_mouseNafld_randomsubset.rds')

annotation.mouseStSt.Subset <- 
  read.csv(file = 'scRNAseq/Liver Cell Atlas Data/annot_mouseStStAll_randomsubset.csv')

annotation.mouseNafld.Subset <- 
  read.csv(file = 'scRNAseq/Liver Cell Atlas Data/annot_mouseNafldAll_randomsubset.csv')

#### StSt 

# if true:
identical(x = annotation.mouseStSt.Subset$cell , 
          y = names(lca.mouseStSt.Subset@active.ident))

identical(x = annotation.mouseNafld.Subset$cell , 
          y = names(lca.mouseNafld.Subset@active.ident))
# then:
Idents(lca.mouseStSt.Subset) <- annotation.mouseStSt.Subset$annot
Idents(lca.mouseNafld.Subset) <- annotation.mouseNafld.Subset$annot

# scRNA-seq QC metric.
lca.mouseStSt.Subset[["percent.mt"]] <- PercentageFeatureSet(lca.mouseStSt.Subset, 
                                                             pattern = "^mt-")

lca.mouseNafld.Subset[["percent.mt"]] <- PercentageFeatureSet(lca.mouseNafld.Subset, 
                                                             pattern = "^mt-")
#jpeg(file="vln_plot.jpeg")
VlnPlot(lca.mouseNafld.Subset, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt" ), 
         ncol = 3,pt.size = F) 

VlnPlot(lca.mouseNafld.Subset, features = "nFeature_RNA") + 
  theme(legend.position = 'none')
VlnPlot(lca.mouseNafld.Subset, features = "nCount_RNA") + 
  theme(legend.position = 'none')
VlnPlot(lca.mouseNafld.Subset, features = "percent.mt", pt.size = F ) + 
  theme(legend.position = 'none')


plot1 <- FeatureScatter(lca.mouseNafld.Subset, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(lca.mouseNafld.Subset, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1
plot2


#### IN THE CLUSTER: subseting, scttransforming, pca and umap were done 

# what do these plots tell us, how should we filter the data? 
# filter
# 
lca.mouseStSt.Subset <- subset(lca.mouseStSt.Subset, 
                               subset = nFeature_RNA > 200 &
                                 nFeature_RNA < 7000 &
                                percent.mt < 10 & 
                                 nCount_RNA < 50000)

lca.mouseNafld.Subset <- subset(lca.mouseNafld.Subset, 
                               subset = nFeature_RNA > 200 &
                                 nFeature_RNA < 3000 &
                                 percent.mt < 7.5 & 
                                 nCount_RNA < 15000)



## scTransform
# single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().\
# lca.mouseStSt.Subset <- SCTransform(lca.mouseStSt.Subset, 
#                                     vars.to.regress = "percent.mt", 
#                                     verbose = FALSE)

##### NOT IN THE CLUSTER
##### load  scTransformed subset (done in cluster) 
lca.mouseStSt.Subset <- readRDS(file = 'scRNAseq/Liver Cell Atlas Data/countTable_mouseStSt_randomsubset_sctransformed.rds')
temp <- is.element(annotation.mouseStSt.Subset$cell, names(lca.mouseStSt.Subset$orig.ident)) 
annotation.mouseStSt.Subset <- annotation.mouseStSt.Subset[temp,]
lca.mouseStSt.Subset@meta.data$sample <- annotation.mouseStSt.Subset$sample
lca.mouseStSt.Subset@meta.data$typeSample <- annotation.mouseStSt.Subset$typeSample
lca.mouseStSt.Subset@meta.data$digest <- annotation.mouseStSt.Subset$digest
lca.mouseStSt.Subset@meta.data$celltype<- annotation.mouseStSt.Subset$annot


lca.mouseNafld.Subset <- readRDS(file = 'scRNAseq/Liver Cell Atlas Data/countTable_mouseNafld_randomsubset_sctransformed.rds')
temp <- is.element(annotation.mouseNafld.Subset$cell, names(lca.mouseNafld.Subset$orig.ident)) 
annotation.mouseNafld.Subset <- annotation.mouseNafld.Subset[temp,]
lca.mouseNafld.Subset@meta.data$sample <- annotation.mouseNafld.Subset$sample
lca.mouseNafld.Subset@meta.data$typeSample <- annotation.mouseNafld.Subset$typeSample
lca.mouseNafld.Subset@meta.data$digest <- annotation.mouseNafld.Subset$digest
lca.mouseNafld.Subset@meta.data$celltype <- annotation.mouseNafld.Subset$annot

### color by different columns of annotation 
#StSt
DimPlot(lca.mouseStSt.Subset, reduction = 'umap',  
        label = T, group.by = 'seurat_clusters' )
DimPlot(lca.mouseStSt.Subset, reduction = 'umap',  
        label = T, group.by = 'digest' )
DimPlot(lca.mouseStSt.Subset, reduction = 'umap',  
        label = T, group.by = 'typeSample' )
DimPlot(lca.mouseStSt.Subset, reduction = 'umap',  
        label = T, group.by = 'sample' )
DimPlot(lca.mouseStSt.Subset, reduction = 'umap',  
        label = T, group.by = 'celltype' )

DimPlot(lca.mouseStSt.Subset, reduction = 'pca', group.by = 'sample')

#Nafld
DimPlot(lca.mouseNafld.Subset, reduction = 'umap', 
        label = T, group.by = 'seurat_clusters' )
DimPlot(lca.mouseNafld.Subset, reduction = 'umap', 
        label = T, group.by = 'digest' )
DimPlot(lca.mouseNafld.Subset, reduction = 'umap', 
        label = T, group.by =  'typeSample' )
DimPlot(lca.mouseNafld.Subset, reduction = 'umap', 
        label = T, group.by = 'sample' )
DimPlot(lca.mouseNafld.Subset, reduction = 'umap',  
        label = T, group.by = 'celltype' )

DimPlot(lca.mouseNafld.Subset, reduction = 'pca', group.by = c('sample','orig.ident'))


# check if cells are clustering 'just because' of their differencens in 
# mt gene expression 
lca.mouseStSt.Subset[["mt.expr.level"]] <- NA
summary(lca.mouseStSt.Subset@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.084   4.348   4.182   6.056   9.998 
# use quantiles to categorize cells in mt.expression
lca.mouseStSt.Subset@meta.data[lca.mouseStSt.Subset@meta.data$percent.mt >= 6.056, "mt.expr.level"] <- 'high'
lca.mouseStSt.Subset@meta.data[lca.mouseStSt.Subset@meta.data$percent.mt >= 2.084 &
                                 lca.mouseStSt.Subset@meta.data$percent.mt <= 6.056, "mt.expr.level"] <- 'medium'
lca.mouseStSt.Subset@meta.data[lca.mouseStSt.Subset@meta.data$percent.mt <= 2.084, "mt.expr.level"] <- 'low'


#  same for Nafld 
lca.mouseNafld.Subset[["mt.expr.level"]] <- NA
summary(lca.mouseNafld.Subset@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.760   3.277   3.328   4.867   7.499 
# use quartiles to categorize cells in mt.expression
lca.mouseNafld.Subset@meta.data[lca.mouseNafld.Subset@meta.data$percent.mt >= 4.867, "mt.expr.level"] <- 'high'
lca.mouseNafld.Subset@meta.data[lca.mouseNafld.Subset@meta.data$percent.mt >= 1.760  &
                                 lca.mouseNafld.Subset@meta.data$percent.mt <= 4.867, "mt.expr.level"] <- 'medium'
lca.mouseNafld.Subset@meta.data[lca.mouseNafld.Subset@meta.data$percent.mt <= 1.760, "mt.expr.level"] <- 'low'


#look at it 
DimPlot(lca.mouseNafld.Subset, reduction = 'umap', group.by = 'mt.expr.level')

# average gene expression 
lca.mouseStSt.Subset.averageExpression <- AverageExpression(object = lca.mouseStSt.Subset,
                                                            assays = 'SCT', group.by = 'celltype')
lca.mouseStSt.Subset.averageExpression <- as.data.frame(lca.mouseStSt.Subset.averageExpression$SCT)
summary(lca.mouseStSt.Subset.averageExpression)


lca.mouseNafld.Subset.averageExpression <- AverageExpression(object = lca.mouseNafld.Subset,
                                                            assays = 'SCT', group.by = 'celltype')
lca.mouseNafld.Subset.averageExpression <- as.data.frame(lca.mouseNafld.Subset.averageExpression$SCT)
summary(lca.mouseNafld.Subset.averageExpression)


# average expression across all cell types except hepatocytes (col 5)
lca.mouseStSt.Subset.averageExpression$global_avg_minus_hep <- 
  rowSums(lca.mouseStSt.Subset.averageExpression[,-8])
lca.mouseNafld.Subset.averageExpression$global_avg_minus_hep <- 
  rowSums(lca.mouseNafld.Subset.averageExpression[,-7])

# Fold Changes

other_cell_types <- unique(annotation.mouseStSt.Subset$annot)
other_cell_types <- other_cell_types[-8]

lca.mouseStSt.Subset.FC.DetectionRate <- FoldChange(
  lca.mouseStSt.Subset,
  group.by = 'celltype',
  ident.1 = 'Hepatocytes',
  ident.2 = other_cell_types,
  assay = 'SCT',
  base = 2,
)

other_cell_types <- unique(annotation.mouseNafld.Subset$annot)
other_cell_types <- other_cell_types[-7]

lca.mouseNafld.Subset.FC.DetectionRate <- FoldChange(
  lca.mouseNafld.Subset,
  group.by = 'celltype',
  ident.1 = 'Hepatocytes',
  ident.2 = other_cell_types,
  assay = 'SCT',
  base = 2,
)


universal_gene_signature <- read.csv(file = 'Universal_Damage_Signature_Kos_16.03.22_Top100.csv',
                   sep = ';')


## this plot I WANT 
ggplot(data= log2(lca.mouseStSt.Subset.averageExpression), 
       aes(x=Hepatocytes,y=global_avg_minus_hep)) + geom_point(alpha = 0.3) + 
  geom_point(data= log2(lca.mouseStSt.Subset.averageExpression[top40_celltype_not_filtered$gene_symbol,]),
                                       aes(x=Hepatocytes,y=global_avg_minus_hep, color = 'Top 40 Genes'), alpha = 0.7)  + xlab('log2(avg. expression per gene in hepatocytes)') +  ylab('log2(avg. expression per gene all other cell types)')  
#StSt_hist_avg_gene_expression_hepatocytes_40top
ggplot(data= log2(lca.mouseStSt.Subset.averageExpression), 
       aes(x=Hepatocytes)) + geom_histogram(binwidth = 0.25) + geom_point(data= log2(lca.mouseStSt.Subset.averageExpression[top40_celltype_not_filtered$gene_symbol,]),
                                                                          aes(x=Hepatocytes,y=0, color = 'Top 40 Genes'), alpha = 0.7) + xlab('log2(avg. expression per gene in hepatocytes)')
#StSt_hist_avg_gene_expression_allothercelltypes_40top
ggplot(data= log2(lca.mouseStSt.Subset.averageExpression), 
       aes(x=global_avg_minus_hep)) +geom_histogram(binwidth = 0.25) + geom_point(data= log2(lca.mouseStSt.Subset.averageExpression[top40_celltype_not_filtered$gene_symbol,]),
                                                                                  aes(x=global_avg_minus_hep,y=0, color = 'Top 40 Genes'), alpha = 0.7) + xlab('log2(avg. expression per gene in all other cell types)')

# for the threshold
summary(log2(lca.mouseStSt.Subset.averageExpression$Hepatocytes))
summary(log2(lca.mouseStSt.Subset.averageExpression$global_avg_minus_hep))

### 2log fold change hepatocytes vs. all other cells 

ggplot(data= lca.mouseStSt.Subset.FoldchangeFuncResults, 
       aes(x=avg_log2FC)) +
  geom_histogram(binwidth = 0.01) + xlim(-3,3) + ylim(0,1000) +
  geom_point(data= lca.mouseStSt.Subset.FoldchangeFuncResults[top40_celltype_not_filtered$gene_symbol,],
             aes(x=avg_log2FC,y=0, color = 'Top 40 Genes'), alpha = 0.6)

# ggplot(data= lca.mouseStSt.Subset.averageExpression, aes(x=Hepatocytes)) +
#   geom_density() + 
#   geom_density(data= lca.mouseStSt.Subset.averageExpression,
#                                        aes(x=global_avg_minus_hep, color = 'Top 40 Genes'), alpha = 0.6) +
#   xlim(-2,10) + ylim(-0.5, 5) + 
#   geom_point(data= lca.mouseStSt.Subset.averageExpression[top40_celltype_not_filtered$gene_symbol,],
#                                            aes(x=Hepatocytes,y=global_avg_minus_hep, color = 'Top 40 Genes'), alpha = 0.6) 
# 

summary(lca.mouseStSt.Subset.FoldchangeFuncResults[top40_celltype_not_filtered$gene_symbol,])
summary(lca.mouseStSt.Subset.FoldchangeFuncResults)

did_not_pass_40_lfc<- lca.mouseStSt.Subset.FoldchangeFuncResults[top40_celltype_not_filtered$gene_symbol,'avg_log2FC'] < -0.047028
did_not_pass_40_lfc<- top40_celltype_not_filtered$gene_symbol[did_not_pass_40_lfc]


ggplot(data= lca.mouseStSt.Subset.averageExpression, 
       aes(x=Hepatocytes)) + geom_histogram(binwidth = 0.1)  + xlim(-0.5,30) + ylim(0,3000) +
  geom_point(data= lca.mouseStSt.Subset.averageExpression[top40_celltype_not_filtered$gene_symbol,],
             aes(x=Hepatocytes,y=0, color = 'Top 40 Genes'), alpha = 0.6)

summary(lca.mouseStSt.Subset.averageExpression$Hepatocytes)

did_not_pass_40_avg<- lca.mouseStSt.Subset.averageExpression[top40_celltype_not_filtered$gene_symbol,'Hepatocytes'] < 0.13538
count(lca.mouseStSt.Subset.averageExpression$Hepatocytes[did_not_pass_40_avg])

did_not_pass_40_avg <- top40_celltype_not_filtered$gene_symbol[did_not_pass_40_avg]

ggplot(data= lca.mouseStSt.Subset.averageExpression, 
       aes(x=Hepatocytes)) + geom_histogram(binwidth = 0.25)  + #xlim(-0.5,30) + ylim(0,3000) +
  geom_point(data= lca.mouseStSt.Subset.averageExpression[top40_celltype_not_filtered$gene_symbol,],
             aes(x=Hepatocytes,y=0, color = 'Top 40 Genes'), alpha = 0.6)  

ggplot(data= log2(lca.mouseStSt.Subset.FoldchangeFuncResults), 
       aes(x=pct.1)) + geom_histogram(binwidth = 0.25)  + #xlim(-0.5,30) + ylim(0,3000) +
  geom_point(data= log2(lca.mouseStSt.Subset.FoldchangeFuncResults[top40_celltype_not_filtered$gene_symbol,]),
             aes(x=pct.1,y=0, color = 'Top 40 Genes'), alpha = 0.6) + xlab('percentage of hepatocytes expressing gene (log2)')

ggplot(data= lca.mouseStSt.Subset.FoldchangeFuncResults, 
       aes(x=avg_log2FC)) + geom_histogram(binwidth = 0.010)  + #xlim(-0.5,30) + ylim(0,3000) +
  geom_point(data= lca.mouseStSt.Subset.FoldchangeFuncResults[top40_celltype_not_filtered$gene_symbol,],
             aes(x=avg_log2FC,y=0, color = 'Top 40 Genes'), alpha = 0.6) + xlab('LFC gene expression in hepatocytes compared to all other cell types')

ggplot(data= log2(lca.mouseStSt.Subset.FoldchangeFuncResults), 
       aes(x=pct.2,y=pct.1)) + geom_point(alpha = 0.3) + 
  geom_point(data= log2(lca.mouseStSt.Subset.FoldchangeFuncResults[top40_celltype_not_filtered$gene_symbol,]),
             aes(x=pct.2,y=pct.1, color = 'Top 40 Genes'), alpha = 0.7)  + xlab('log2(avg. expression per gene in hepatocytes)') +  ylab('log2(avg. expression per gene all other cell types)')  

