##############################################################
### calculating damage score, using a gene set and AUCell ####
########################################################## ####

### make a collection of gene sets  to use with AUCell
library(GSEABase)
# genesUP and genesDOWN are charachter vectors containing 
# the names of genes UP and DOWN regulated upon damage 
genes <- read.csv('Universal_Damage_Signature_Kos_18.03.22_conservative_list.csv', sep = ';')
genesUP <- genes[genes$direction_foldchange == 1, 'gene_symbol']
genesDOWN <- genes[genes$direction_foldchange == -1, 'gene_symbol']

geneSets <-  GSEABase::GeneSetCollection( list( GeneSet( genesUP , 
                                               setName= "UP" ) , 
                                      GeneSet( genesDOWN , 
                                               setName= "DOWN" )) )


##  LOAD DATASET to apply damage score to:

# 1) single cell RNAseq: 

library(Seurat)
# Getting scRNAseq data set (seurat object) ready
load("/Users/paulaungera/Desktop/Project Module Beyer /scRNAseq/GSE148859_SingleCellExprement.Rdata")
sce.combined[["percent.mt"]] <- PercentageFeatureSet(sce.combined, pattern = "mt-")
VlnPlot(sce.combined, features = "percent.mt")
#no filtering by mt DNA because all cells have SO MUCH of it, I end up with too little cells 

sce.combined <- NormalizeData(sce.combined)
#check out cell identities 
Idents(sce.combined)
DimPlot(sce.combined, reduction = "umap")

subset_hepatocytes <- subset(x = sce.combined, 
                       idents = c("Hepatocytes 1",
                                  "Hepatocytes 2",
                                  "Hepatocytes 3",
                                  "Hepatocytes 4"))

exprMatrices <- GetAssayData(object = subset_hepatocytes , slot = 'counts')
exprMatrices <- as.matrix(exprMatrices)
saveRDS(exprMatrices, file = 'scRNAseq/GSE148859_SingleCellExperiment_exprMatrices.rds')



# or
## 2) FOR BULK
# one of 'my' bulk sets for testing 
# exprMatrices <- read.csv('GSE137449/GSE137449_count_table', sep=' ')
# rownames(exprMatrices) <- exprMatrices$gene_symbol
# exprMatrices <- exprMatrices[,-1]
# exprMatrices <- as.matrix(exprMatrices)


require(AUCell)

###  1. Build gene-expression rankings for each cell 
cellRanks <- AUCell_buildRankings(  exprMatrices, nCores=1, plotStats=F)

###  2. Calculate enrichment for the gene signatures (AUC)
# gene sets should include UP and DOWN regulated genes seperately
cells_AUC <- AUCell_calcAUC( geneSets, cellRanks , 
                             # aucMaxRank parameter is a tricky one,
                             # for single cell 0.05 is a reasonable 
                             # but for bulk a bigger, e.g. 0.2 can be used
                             aucMaxRank = ceiling(0.05 * nrow(cellRanks))
                             )

cells_AUC <- getAUC( cells_AUC) #this will extract results 

### calculate the damage score
# indeces specify elements of geneSets that contain UP and DOWN regulated genes, corespondingly
index1 <- 1
index2 <- 2
Dscore= cells_AUC[index1,]-cells_AUC[index2,] 


##### PLOT RESULTS ######
library(ggplot2)
library(ggpubr)

# 1) single cell RNAseq:

#### Example scRNAseq GSE148859: plot Dscore

temp <-data.frame('samples' =names(Dscore), 'Dscore' = as.vector(Dscore))

# here you sort out how you want to group the cell in which boxplots
identity <- subset_hepatocytes@meta.data$orig.ident
# results in 1-254 WT
identity <- sub('_1', '', identity)
identity <- sub('_2', '', identity)
identity <- sub('_3', '', identity)
temp$genotype <- NA 
temp$genotype <- identity
temp$genotype <- gsub(pattern = 'HOHO', 'Tak1_HEP+D138', temp$genotype)
temp$genotype <- gsub(pattern = 'HO', 'Tak1_HEP', temp$genotype)
unique(temp$genotype)

#plot
p <- ggboxplot(temp, x = 'genotype', y = "Dscore",
               color = 'genotype', palette = "jco",
               add = 'jitter')

#  Add p-value (Wilcoxon)
p + stat_compare_means()



# 2) bulk RNAseq: 
annotation <- read.csv('GSE137449/SraRunTable.txt')
temp <-data.frame('samples' =names(Dscore), 'Dscore' = as.vector(Dscore))
temp$treatment <- annotation[annotation$Genotype == 'wild type (WT)', 'diet']
temp[temp$treatment != 'normal chow', 'treatment'] <- 'CDAHFD'

p <- ggboxplot(temp, x = 'treatment', y = "Dscore",
               color = "treatment", palette = "jco",
               add = "jitter")

#  Add p-value (Wilcoxon)
p + stat_compare_means()


