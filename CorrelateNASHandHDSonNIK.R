library(AUCell)
library(BiocGenerics)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(plyr)
library(biomaRt)
library(Seurat)
library(ggpubr)
library(cowplot)
library(magrittr)
library(tidyverse)
library(dplyr)
library(psych)
library(tidyr)
library(ggdendro)
library(scales)
library(clusterProfiler)

#### load functions written by us
# generate damage signature
{
  source("~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
  # calculate damage score
  source("~/Repositories/cell-damage-score/AUCell_script.r")
}

pathNIK <- '~/Desktop/HepatocyteDamageScore/scRNAseq/Thomas_Wunderlich_snRNAseq/nik_samples_annotated.rds' 

outputPath <- '~/Desktop/HepatocyteDamageScore/Results/CorrelateNashHDS/'

NIKObject <- readRDS(pathNIK) %>%
  subset(., subset = class_annotation == "Hepatocyte") 

annotation <- NIKObject@meta.data

## Is this good?  
exprMatrix <-  GetAssayData(object = NIKObject , slot = 'counts')  
exprMatrix <- exprMatrix[ rowSums( exprMatrix > 0 ) > ncol(exprMatrix)*0.005 , ]

backgroundForGSEA <- rownames(exprMatrix)

HUDS <- read.csv(file = 
                   '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')

# needs to be replaced by final version
geneFilter <- read.csv(file ='~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')

# 1. Calculate HDS for all hepatocytes
#devide cells into 5 groups
chunks <- seq(dim(exprMatrix)[2]/5 , dim(exprMatrix)[2], dim(exprMatrix)[2]/5 )
i <- 0
NikHDS <- DS_calc.func(exprMatrices = exprMatrix[,1:chunks[1]],
                       DSignature = HUDS, 
                       useOrder = "mean_rank",
                       ntop = 42)
# for( i in 1:4){
#   NikHDS <-c( NikHDS, DS_calc.func(exprMatrices = exprMatrix[,(chunks[i]+1):chunks[i+1]],
#                                    DSignature = HUDS, 
#                                    useOrder = "mean_rank",
#                                    ntop = 42))
# }

# reformat HDS + annotations
NikHDS <- data.frame( 'HDS' = NikHDS, "cellID" = names(NikHDS))
NikHDS$condition <- gsub(pattern = ".*_", replacement = "", NikHDS$cellID)
NikHDS$genotype <- gsub(pattern = ".*_NIK", replacement = "", NikHDS$cellID) %>% 
  gsub(pattern = "_.*", replacement = "", NikHDS$cellID)
NikHDS$cellID <- gsub(pattern = "_.*", replacement = "", NikHDS$cellID)
rownames(NikHDS) <- 1:length(NikHDS$HDS)


## don't know if data that I have is transformed?? 
corList <- list()
conditionsNames <- c("NIKKO_Nash","NIKFL_Nash", "NIKKO_Cdaa","NIKFL_Cdaa", "NIKKO_Hcc", 
                      "NIKFL_Hcc")
for (conditionPattern in conditionsNames){
  
  subsetIndices <- grep(x = colnames(exprMatrix), 
                        pattern = conditionPattern )
  
  corList[[conditionPattern]] <- 
    psych::corr.test(x = as.matrix(t(exprMatrix[, subsetIndices])),
                     y = NikHDS$HDS[subsetIndices],
                     method = 'spearman')
  
  
}

# saveRDS(corList, file = '~/Desktop/HepatocyteDamageScore/Results/CorrelationHDStoExprNIKPerCondition.rds')

corList <- readRDS('~/Desktop/HepatocyteDamageScore/Results/CorrelationHDStoExprNIKPerCondition.rds')

# heatmap genes correlation expr and HDS vs. condition (not what Tim wanted, 
# he wanted each column of heatmap to be a sample)
corDataFrame <- data.frame(corList[[1]]["r"], 
                           corList[[2]]["r"], 
                           corList[[3]]["r"],
                           corList[[4]]["r"], 
                           corList[[5]]["r"],
                           corList[[6]]["r"]
                             )
colnames(corDataFrame) <- names(corList)
corDataFrame$genes <- row.names(corDataFrame)


sigGenes <- list()
allGenes <- rownames(corList[[1]][["p"]])
sigGenes <- lapply(1:6, function(ii){
  
  sigGenes[[ii]] <- allGenes[which(corList[[ii]][["p"]] <= 0.01 &
                                     abs(corList[[ii]][["r"]]) >= 0.2)]
  
}
)

sigGenesUnion <- Reduce(f = union, x = sigGenes)

# replace NAs with zeroes
toPlot <- corDataFrame[match(sigGenesUnion, corDataFrame$genes),-7] %>% replace(is.na(.), 0) 
rownames(toPlot) <- sigGenesUnion
# dataframe with correlation between HDS and gene expression
# not significantly genes correlated genes filtered out (union of significant)
# genes per model
# toPlot <- data.frame("NIKFL_Nash" = toPlot$NIKFL_Nash,
#                      "NIKKO_Nash" = toPlot$NIKKO_Nash,
#                      "NIKFL_CDAA" = toPlot$NIKFL_Cdaa,
#                      "NIKKO_CDAA" = toPlot$NIKKO_Cdaa,
#                      "NIKFL_HCC" = toPlot$NIKFL_Hcc,
#                      "NIKKO_HCC" = toPlot$NIKKO_Hcc)

HeatMapClustered <- 
  pheatmap::pheatmap(toPlot,
                     show_rownames = FALSE,
                     cluster_cols = FALSE,
                     color = colorRampPalette(c("blue","white", "red"))( 50 ),
                     breaks = seq(-0.7, 0.7, 1.4/50),
                     clustering_method = 'ward.D2'
                     ) 

png(paste0(outputPath,'HeatmapClusterSigCorrGenesHDS.png'), 
    units = 'px', height = 700, width =  600, )
HeatMapClustered
dev.off()

## subtract correlation values from FL and KO to highlight differences and plot
# heatmap with clustering

deltaGeno <- data.frame("deltaNash" = toPlot[,1] - toPlot[,2],
                        "deltaCdaa" = toPlot[,3] - toPlot[,4],
                        "deltaHcc" = toPlot[,5] - toPlot[,6])

## so here if 
HeatMapClusteredGeno <- 
  pheatmap::pheatmap(deltaGeno,
                     show_rownames = FALSE,
                     cluster_cols = FALSE,
                     color = colorRampPalette(c("blue","white", "red"))( 50 ),
                     breaks = seq(-0.7, 0.7, 1.4/50),
                     clustering_method = 'ward.D2'
  ) 

## subtract correlation values from HCC to the other models to highlight 
# differences and plot
# heatmap with clustering --> end with four cols

deltaHcc <- data.frame("deltaNashFL" = toPlot[,1] - toPlot[,5],
                       "deltaNashKO" = toPlot[,2] - toPlot[,6],
                        "deltaCdaaFL" = toPlot[,3] - toPlot[,5],
                       "deltaCdaaKO" = toPlot[,4] - toPlot[,6]
                       )

## so here if 
summary(deltaHcc)
HeatMapClusteredDeltaHcc <- 
  pheatmap::pheatmap(deltaHcc,
                     show_rownames = FALSE,
                     cluster_cols = FALSE,
                     color = colorRampPalette(c("blue","white", "red"))( 50 ),
                     breaks = seq(-0.4, 0.4, 0.8/50),
                     clustering_method = 'ward.D2'
  ) 




# if you have two samples per group you can do differential 
# testing per group and use a linear model to see which genes 
# contribute to separate diet and genotype of an interaction of both 

## Tim's clustering:


# 1. cluster your data using with the method of your choice
# dist() calculated euclidean distance matrix, then hclust clusters 
# the genes according to their distance in distance matrix
clust.obj <- hclust( dist(toPlot), method="ward.D2")
# with clust.obj you can make trees, heatmaps, etc. 
# 2. cut tree to get N clusters (k), or at a certain height (h)
N = 8
cclust <- cutree( clust.obj , k = N) 
cclust <- as.factor( cclust )
# 3. assign colors to your clusters
levels(cclust) <- hue_pal()(N)
# 4. plot a heatmap with colored clusters
pdf(file = paste0(outputPath, 'HeatmapNikGeneExprHDSCorrelation.pdf'))

heatmap <- gplots::heatmap.2( t(as.matrix(toPlot)),
                   labCol =""  ,
                   trace = "none" , 
                   Rowv = FALSE ,
                   dendrogram = 'column',
                   col = colorRampPalette(c("blue","white", "red"))( 50 ), 
                   # here chose whatever colors you want
                   hclustfun = function(x) hclust(x, 
                                                  method="ward.D2") , 
                   # should be the same clustering  method as in step 1!
                     ColSideColors= as.character( cclust ) , 
                   RowSideColors=rep(c("orange","brown", "black"), each=2),
                   margin=c(9, 9),
                   # this will add a cluster annotation bar
                   main = "HDS and gene expression\n correlation")



dev.off()
unique(cclust)

# colors of clusters of interest
# have to check if this is right order or the other way around 
clustersInterest <- c("#7CAE00",'#C77CFF','#FF61CC')

ClustersDataFrame <- data.frame('gene' = names(cclust),
                                'cluster' = cclust,
                                row.names = 1:length(cclust))
write.csv(ClustersDataFrame[ClustersDataFrame$cluster == clustersInterest[1], 
                            'gene'], 
          file = paste0(outputPath, 'OliveGreenClusterGenesHDSNIK.csv'),
          quote = FALSE,
          row.names = FALSE)

write.csv(ClustersDataFrame[ClustersDataFrame$cluster == clustersInterest[2], 
                            'gene'], 
          file = paste0(outputPath, 'PurpleClusterGenesHDSNIK.csv'),
          quote = FALSE,
          row.names = FALSE)

write.csv(ClustersDataFrame[ClustersDataFrame$cluster == clustersInterest[3], 
                            'gene'], 
          file = paste0(outputPath, 'PinkClusterGenesHDSNIK.csv'),
          quote = FALSE,
          row.names = FALSE)

# background

write.csv(backgroundForGSEA, 
          file = paste0(outputPath, 'CorrHDSNIK_BackgroundGenesForGSEA.csv'),
          quote = FALSE,
          row.names = FALSE)

write.csv(ClustersDataFrame,
          file = paste0(outputPath, 'AllClustersGenesHDSNIK.csv'),
          quote = FALSE,
          row.names = FALSE)

