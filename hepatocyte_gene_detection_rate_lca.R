library(Seurat)
library(sctransform)
library(ggplot2)
## merge

temp_stst <- readRDS(file = '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/countTable_mouseStSt_Hepatocytes_nucSeq.rds')
annotation_temp <- read.csv(file ='~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/annot_mouseStStAll.csv')
annotation_temp <- annotation_temp[annotation_temp$annot == 'Hepatocytes', ]
annotation_temp<- annotation_temp[is.element(annotation_temp$cell, names(temp_stst$orig.ident)),] 
temp_stst@meta.data$sample <- annotation_temp$sample
temp_stst@meta.data$digest <- annotation_temp$digest
temp_stst@meta.data$typeSample <- annotation_temp$typeSample
temp_stst@meta.data$celltype<- annotation_temp$annot

temp_nafld <- readRDS(file = '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/countTable_mouseNafld_Hepatocytes_nucSeq.rds')
annotation_temp <- read.csv(file ='~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/annot_mouseNafldAll.csv')
annotation_temp <- annotation_temp[annotation_temp$annot == 'Hepatocytes', ]
annotation_temp<- annotation_temp[is.element(annotation_temp$cell, names(temp_nafld$orig.ident)),] 
temp_nafld@meta.data$sample <- annotation_temp$sample
temp_nafld@meta.data$digest <- annotation_temp$digest
temp_nafld@meta.data$typeSample <- annotation_temp$typeSample
temp_nafld@meta.data$celltype<- annotation_temp$annot


merged.mouse.lca.hep.nucSeq <-
  merge(x = temp_stst,
        y = temp_nafld, 
        add.cell.ids = c("StSt", "Nafld"), 
        project = "LiverCellAtlas")

remove(temp_nafld,temp_stst, annotation_temp)

saveRDS(merged.mouse.lca.hep.nucSeq, file = 'countTable_mouseMerged_Hepatocytes_nucSeq.rds')

#Quality control 

merged.mouse.lca.hep.nucSeq <- readRDS('~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/countTable_mouseMerged_Hepatocytes_nucSeq.rds')

merged.mouse.lca.hep.nucSeq[["percent.mt"]] <- PercentageFeatureSet(merged.mouse.lca.hep.nucSeq, 
                                                                    pattern = "^mt-")

# check dimension before and after next step
dim(merged.mouse.lca.hep.nucSeq)
# 26269 49521

merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, subset = typeSample == 'nucSeq' )
# after: 26269 42927


VlnPlot(merged.mouse.lca.hep.nucSeq, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt" ), 
        ncol = 3,pt.size = F) 


VlnPlot(merged.mouse.lca.hep.nucSeq,
        features = c("percent.mt" ), pt.size = F, group.by = 'typeSample' ) 


plot1 <- FeatureScatter(merged.mouse.lca.hep.nucSeq, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(merged.mouse.lca.hep.nucSeq, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1
plot2
## maybe cluster not necessary because small
## Choose filtering parameters
# nFeature_RNA > ... & nFeature_RNA < ... & percent.mt < ... & nCount_RNA < ...
## filter in the cluster and sct transform and find neighbors and clusters
# also dimension reduction
# reload data 

#only single nuclei data



merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, 
                               subset = nFeature_RNA > 200 &
                                 nFeature_RNA < 6000 &
                                 percent.mt <= 0 & 
                                 nCount_RNA < 30000 )

# this numbers down here are not accurate, they probably belong to less stringent
# filtering
#liver_cell_atlas_Nafld  liver_cell_atlas_stst 
#20775                  10829 

# this is what I get:
# 26269  3291

# [1] 26269  2567 numbers now

merged.mouse.lca.hep.nucSeq <- SCTransform(merged.mouse.lca.hep.nucSeq,
                                     verbose = FALSE)


# aprox. number of dimension adter elbow plot
merged.mouse.lca.hep.nucSeq <- RunPCA(merged.mouse.lca.hep.nucSeq, 
                                      verbose = FALSE)
ElbowPlot(merged.mouse.lca.hep.nucSeq)

merged.mouse.lca.hep.nucSeq <- FindNeighbors(merged.mouse.lca.hep.nucSeq, 
                                             dims = 1:10)
merged.mouse.lca.hep.nucSeq <- FindClusters(merged.mouse.lca.hep.nucSeq, 
                                            resolution = 0.5)
merged.mouse.lca.hep.nucSeq  <- RunUMAP(merged.mouse.lca.hep.nucSeq , 
                                        dims = 1:10, verbose = FALSE)

# plot dimension redution control 
# 0. by treament
DimPlot(merged.mouse.lca.hep.nucSeq, reduction = 'umap',
        label = T, group.by = 'orig.ident')
DimPlot(merged.mouse.lca.hep.nucSeq, reduction = 'pca',
        group.by = 'orig.ident')

# 1. by cluster
DimPlot(merged.mouse.lca.hep.nucSeq, reduction = 'umap',  
        label = T, group.by = 'seurat_clusters')
#2. by samples
DimPlot(merged.mouse.lca.hep.nucSeq, reduction = 'umap',  
        label = T, group.by = 'sample')
# compare to sample info

sample_info <- read.table(file ='~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/GSE192740_sampleInfo_scRNAseq.tsv',
                          header = T, sep = '\t')

# I don't know what I wanted to do here :D
sample_info <- sample_info[sample_info$molecule == 'Nuclear RNA',]
table(sample_info$characteristics..strain)
temp <- unlist(merged.mouse.lca.hep.nucSeq$sample)
merged.mouse.lca.hep.nucSeq$strain <- merged.mouse.lca.hep.nucSeq$sample

temp <- sapply(temp, function(i) i <- sample_info["characteristics..shortFileName" == i,"characteristics..strain" ] )


## calculate parameters for gene filtering

# Detection Rate in hepatocytes (LFC between single nucleus and single cell)
foldchange <- FoldChange(merged.mouse.lca.hep.nucSeq,
                         ident.1 = 'nucSeq', 
                         group.by = 'typeSample',
                         assay = "SCT")

detection_rate_df <- data.frame('Genes' = row.names(foldchange), 
                             'detection_rate' = foldchange[,2],
                             'log2_dec_rate' = log2(foldchange[,2]+0.0001))

ggplot(data= detection_rate_df, aes(x=detection_rate)) + geom_histogram(binwidth = 0.005)

summary(detection_rate_df$log2_dec_rate, na.rm=T)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-13.287712  -8.333516  -6.493297  -6.516940  -4.714065   0.000144 

quantile(detection_rate_df$log2_dec_rate)
#0%           25%           50%           75%          100% 
#-1.328771e+01 -8.333516e+00 -6.493297e+00 -4.714065e+00  1.442623e-04 
quantile(detection_rate_df$detection_rate)
#0%   25%   50%   75%  100% 
#0.000 0.003 0.011 0.038 1.000

plot <- ggplot(data= detection_rate_df, aes(x=log2_dec_rate)) +
  geom_histogram(binwidth = 0.10) + geom_vline(xintercept = c(-8.333516,-6.493297,-4.714065), colour = 'red') +
  xlab('log2(detection rate of gene + 0.0001)') + ylab('number of genes') + 
  geom_text(aes(x=-8.333516, label="1. Quartile", y=1100), colour="red", angle=90, vjust = -1) +
  geom_text(aes(x=-6.493297, label="2. Quartile", y=1100), colour="red", angle=90, vjust = -1) +
  geom_text(aes(x=-4.714065, label="3. Quartile", y=1100), colour="red", angle=90, vjust = -1)

detection_rate_df_selected <- detection_rate_df[detection_rate_df$log2_dec_rate >  -5,]

detection_rate_df_selected <- detection_rate_df_selected[order(detection_rate_df_selected$log2_dec_rate, decreasing = T),]
rownames(detection_rate_df_selected) <- 1:length(detection_rate_df_selected$Genes)

# 
write.csv(detection_rate_df_selected, file = '~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/p_nucseq_genes_filtered_genelist_minus5threshold_log2TESTtoCOMPARE.csv')

# write_csv(detection_rate_df_selected, file = 'p_nucseq_genes_filtered_genelist_minus5threshold_log2.csv')

# test to compare: 
geneSetNow <- detection_rate_df_selected$Genes
genesSetThen <- read.csv(file = '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/p_nucseq_genes_filtered_genelist_minus5threshold_log2.csv')
genesSetThen <- genesSetThen$Genes

length(intersect(geneSetNow, genesSetThen))
length(genesSetThen)
