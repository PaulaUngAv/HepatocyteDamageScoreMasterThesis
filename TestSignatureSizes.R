# ============================================#
############# Test Signature Sizes ############ 
# ============================================#

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

#### load functions written by us
# generate damage signature
{
source("~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
# calculate damage score
source("~/Repositories/cell-damage-score/AUCell_script.r")
}

#### Load Count Data
{
pathToData <- '~/Desktop/HepatocyteDamageScore/'
bulkStudies <- c('GSE99010/GSE99010_CCl4_deg.csv',
             'GSE99010/GSE99010_WesternDiet_deg.csv',
             'GSE97234/GSE97234_AH_deg.csv',
             'GSE97234/GSE97234_ASH_deg.csv',
             'GSE83240/GSE83240_DEN_deg.csv',
             'GSE153580/GSE153580_STZ_deg.csv',
             'GSE148849/GSE148849_fastfooddiet_deg.csv',
             'GSE138419/GSE138419_AMLN_deg.csv',
             'GSE137449/GSE137449_CDAHFD_diet_deg.csv',
             'GSE135050/GSE135050_hfcfdiet_deg.csv',
             'GSE132040/GSE132040_young(6_9_12)vs_old(18_21_24)_deg.csv',
             'GSE123894/GSE123894_Fructose_only_DBA2Jstrain_deg.csv',
             'GSE119953/GSE119953_HFD_deg.csv',
             'GSE119953/GSE119953_DDC_deg.csv',
             'GSE119441/GSE119441_HFD_deg.csv',
             'GSE119441/GSE119441_PFOA_deg.csv',
             'GSE114261/GSE114261_STAM_deg.csv',
             'GSE111828/GSE111828_Acetaminophen_deg.csv'
)

data <- readRDS(
  '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv')
annotationPath <- 
  '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv'
annotation <- read.csv(annotationPath, row.names = 1) 
sampleInfo <- read.csv('~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data//GSE192740_sampleInfo_modified.csv', sep = ";")
# 
# scExprMatrices <- c('~/Desktop/HepatocyteDamageScore/scRNAseq/GSE148859/GSE148859HepatocytesExprMatrix.rds',
#                 '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv',
#                 '~/Desktop/HepatocyteDamageScore/scRNAseq/scRNAseq/Thomas_Wunderlich_snRNAseq/nik_samples_annotated.rds')
# sc_annotations <- c('scRNAseq/GSE148859/sc_annotations_GSE148859_hepatocytes.csv',
#                     '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv')

# loads list of all DEG tables
DElist <- lapply( seq( bulkStudies ), function(ii){
  read.csv(paste0(pathToData,"AllBulkStudiesDEGAnalysis/", bulkStudies[[ii]]))
})
names(DElist) <- basename(bulkStudies)

# Hepatocyte Universal Damage Signature
HUDS <- read.csv(file = 
                   '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')

# needs to be replaced by final version
geneFilter <- read.csv(file ='~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')

}

#### Test Size on Liver Cell Atlas Hepatocytes #### 
# only single nuclues RNA seq hepatocytes
### load expression matrices
# exprMatrix <- readRDS(scExprMatrices[2])
# annotation <- read.csv(sc_annotations[2], row.names = 1) 
# set.seed(1)
# keep only NAFLD samples 24 weeks old
annotation$diet <- NA
annotation[annotation$Sample == "CS197" | 
             annotation$Sample ==  "CS196" |
             annotation$Sample == "CS193" |
             annotation$Sample == "CS192" |
             annotation$Sample == "CS191" |
             annotation$Sample == "CS185"  |
             annotation$Sample == "CS183" , 'diet'] <- 'Western Diet'
annotation$diet[is.na(annotation$diet)] <- 'Standard Diet'

toKeep <- annotation$cell_id[annotation$Sample == "CS193" | 
                               annotation$Sample ==  "CS192" |
                               annotation$Sample == "CS191" |
                               annotation$Sample == "CS185" |
                               annotation$Sample == "CS183" |
                               annotation$Sample == "CS190" |
                               annotation$Sample == "CS189" |
                               annotation$Sample == "CS188" |
                               annotation$Sample == "CS184" |
                               annotation$Sample == "CS182"] 
  
  
  #random sampling
   # c(sample(annotation$cell_id[annotation$condition == 'StSt'],
   #          size = 2500, replace = FALSE),
   #   sample(annotation$cell_id[annotation$condition == 'Nafld'], 
   #          size = 2500, replace = FALSE)
   # )
 exprMatrix <- data[, toKeep]
# 

{
  
  size_vec <- c(5, 10, 15, 20, 30 , 42, 50, 100, 200, 300, 500, nrow(HUDS))
  
  ### analyse effect of a set size on individual studies
  
  {
    # iterate over set sizes and calculate the disease score
    
    DSsize <- lapply( seq(size_vec) , 
                      
                       function(ii)
                       {
                         print(size_vec[ii])
                         ntop <- size_vec[ii]
                         DS_calc.func( exprMatrices = exprMatrix,
                                       DSignature = HUDS[1:ntop,]
                         )
                       } )
    
    
    # return DS
    return( DSsize) 
  }
 
  
  # save results in RDS file containing list of names vectors 
  
  # saveRDS(DSsize, 
   #       '~/Desktop/HepatocyteDamageScore/Results/SignatureSizeTest/testSizesmouseMergedHepatocytesnucSeqQCfitleredSCT.rds')
  
  #DSsize <- readRDS('~/Desktop/HepatocyteDamageScore/Results/SignatureSizeTest/testDSsignatureSizesLiverCellAtlas.rds')
  
  ### Prepare data for plotting
  
  DSsizeListDataFrames <- 
    lapply( seq(DSsize) , function(ii){
      datt <- cbind.data.frame( score = DSsize[[ii]] ,
                              setSize = as.factor(size_vec[ii]), 
                              condition = as.factor(
                                annotation$diet[match(names(DSsize[[ii]]),
                                                      annotation$cell_id)])
                              )
    datt$condition <- relevel( datt$condition , ref = "Standard Diet")
    
    return(datt)
  })
  
####  Plotting  
  # about the plotting functions:
  # The function mean_sdl is used for adding mean and standard deviation. 
  # It computes the mean plus or minus a constant times the standard 
  # deviation. In the R code above, the constant is
  # specified using the argument mult (mult = 1). 
  # By default mult = 2. The mean +/- SD can be added as a crossbar or
  # a pointrange.
{
  DSsizeListViolin <- lapply( 
    seq(DSsizeListDataFrames) , function(ii){
      ppg <- ggplot( DSsizeListDataFrames[[ii]] , 
                     aes( x=condition , y=score,
                          color=condition)) + 
        geom_violin(trim = FALSE) + #theme_bw() + 
        scale_color_colorblind() +
        ggtitle(paste("Top", DSsizeListDataFrames[[ii]]$setSize[1], "Genes")) +
        theme( text = element_text(size=12),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.x = element_blank(),
               plot.title = element_text(size=12),
               legend.position = 'none') +
        guides(colour = guide_legend(override.aes = aes(label = "")))+
        labs(y = "HDS") + 
        stat_summary(fun.data = "mean_sdl",
                     fun.args = list(mult = 1),
                     geom = "pointrange"
                     ) 
      return( ppg )
    }
  )
}

SizePannelFigure <- (DSsizeListViolin[[1]] | DSsizeListViolin[[2]] | DSsizeListViolin[[3]] | DSsizeListViolin[[4]]) / (DSsizeListViolin[[5]] | DSsizeListViolin[[6]] | DSsizeListViolin[[7]] | DSsizeListViolin[[8]]) / (DSsizeListViolin[[9]] | DSsizeListViolin[[10]] | DSsizeListViolin[[11]] | DSsizeListViolin[[12]] + theme(legend.position = 'right'))

ggsave(filename = 'SizeTestLiverCellAtlasPanelFigureOnly24WeeksNAFLDCohort.png', 
       units = 'px',
       path = '~/Desktop/HepatocyteDamageScore/Results/SignatureSizeTest/',
       width = 2600, height = 2800)
}
