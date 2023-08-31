#### ==================#####
#### Randomization Test #### 
#### ==================#####

#### Load packages, data and functions ####
{
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
  
  
  # load functions written by us
  # generate damage signature
  source("~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
  # calculate damage score
  source("~/Repositories/cell-damage-score/AUCell_script.r")
  
  # paths to DEG tables of bulk studies
  
  path <- '~/Desktop/HepatocyteDamageScore/'
  studies <- c('GSE99010/GSE99010_CCl4_deg.csv',
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
  # 
  # scExprMatrices <- 
  #   c('~/Desktop/HepatocyteDamageScore/scRNAseq/GSE148859/GSE148859HepatocytesExprMatrix.rds',
  #     '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv', 
  #     '~/Desktop/HepatocyteDamageScore/scRNAseq/scRNAseq/Thomas_Wunderlich_snRNAseq/nik_samples_annotated.rds')
  # 
  #  sc_annotations <- 
  #    c('scRNAseq/GSE148859/sc_annotations_GSE148859_hepatocytes.csv', 
  #      '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv')
  # 
  data <- readRDS(
    '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv')
  annotationPath <- 
    '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv'
  annotation <- read.csv(annotationPath, row.names = 1) 
  sampleInfo <- read.csv('~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data//GSE192740_sampleInfo_modified.csv', sep = ";")
  # 
  # loads list of all DEG tables
  DElist <- lapply( seq(studies), function(ii){
    read.csv(paste0(path,"AllBulkStudiesDEGAnalysis/",studies[[ii]]))}
    )
  names(DElist) <- basename(studies)
  
  # !! need to be replaced by final version
  HUDS <- 
    read.csv(file = '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')
  
  # !!needs to be replaced by final version
  geneFilter <- 
    read.csv(file ='~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')
  geneFilter <- geneFilter$Genes
  
  # Load Liver Cell Atlas Hepatocytes expression matrix for testing #  
  # only single nuclues RNA seq hepatocytes
  # exprMatrix <- readRDS(scExprMatrices[2])
  # annotation <- read.csv(sc_annotations[2], row.names = 1) 
  # sampling 2000 random nuclei per condition 
  # set.seed(2)
  # use only cells 24 weeks NAFLD cohort
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
    #          size = 1000, replace = FALSE),
    #   sample(annotation$cell_id[annotation$condition == 'Nafld'], 
    #          size = 1000, replace = FALSE)
    # )
  exprMatrix <- data[, toKeep]
}

#### Randomization Step ####   
# Randomize a percentage of p-values and LFCs in DEGs and select genes 
# for HUDS, then calculate HDS and compare conditions.
perturbationTestFunction <- function(pert.percent = seq(10, 100, 10),
                                     ntop = 42, exprDat, DEdat, k = 50, 
                                     hepGenes){
  
   
  perturb_test <-lapply( seq(pert.percent) , function(ii) {

    
    print(paste("percent perturbed", pert.percent[ii] ))
    print(paste("these many studies:", length(DEdat)*( pert.percent[ii]/100)))
          # subset data 50 times to avoid data set bias (k times)

          HDS_pert.lvl <- Reduce( cbind.data.frame , lapply( 1:k, function(kk)
          {
            
            print(paste("randomisation round", kk ))
            
            datt <- DEdat
            # choose dataset to randomize
            whichR <- sample( 1:length(datt), length(datt)*( pert.percent[ii]/100) )
            
            # randomize p-value and lfc of chosen datasets
            print(whichR)
            datt[whichR]  <- lapply( datt[whichR]  , function(XX){
              XX$log2FoldChange <- sample( XX$log2FoldChange)
              XX$pvalue <- sample( XX$pvalue)
              return(XX) })
            
            HDS <- damage_signature.func( datt, 
                                          geneFilter = hepGenes)
            
            # calculate disease score for one dataset, sc Liver Cell Atlas, 2000 cells per condition
            XXX <- DS_calc.func( exprMatrices = exprDat, 
                                 DSignature = HDS, ntop = ntop)
            
            
            return(XXX)
          } )) 
          
          # prepare for plotting with ggplot
          # XX <- reshape::melt.list( data=HDS_pert.lvl )
          # XX$random.percent <- as.factor(pert.percent[ii])
          # return(XX)
          HDS_pert.lvl$cell <- rownames(HDS_pert.lvl)
          XX <- reshape::melt( data=HDS_pert.lvl , id.vars = "cell")
          XX$random.percent <- as.factor(pert.percent[ii])
          return(XX)
  })
    
  
    # returns a list of dataframes  
   return(perturb_test) 
    
}

testResult <- perturbationTestFunction(pert.percent = seq(0,100, 10),
                         ntop = 42, exprDat = exprMatrix, 
                         DEdat = DElist, k = 50, hepGenes = geneFilter)


saveRDS(testResult, file = '~/Desktop/HepatocyteDamageScore/Results/RandomizationTest/RandomizationLCAmergedNuq24WeeksNAFLDCohort.rds')

# testResult <- readRDS('~/Desktop/HepatocyteDamageScore/Results/RandomizationTest/RandomizationLCAmergedNuq24WeeksNAFLDCohort.rds')

# Format Results for better plotting
toPlot <- Reduce( rbind, lapply( testResult , function(X) {
  X$value <- scale(X$value, scale = FALSE)
  return(X)
}))

# toPlot$condition <- gsub(pattern = "_.*", replacement = "", toPlot$cell)
# toPlot$cell_id <- gsub(pattern = ".*_", replacement = "", toPlot$cell)
 toPlot$HDS <- toPlot$value
 toPlot <- toPlot[,c(-2, -3)]

#idea 
toPlot$condition <- as.factor(annotation$diet[match(toPlot$cell,
                        annotation$cell_id)])

toPlot$condition <- factor(toPlot$condition,
                           ordered = TRUE,
                           levels = c("Standard Diet", "Western Diet"))

# plotting 
ggp <- ggplot( toPlot , aes( x=random.percent , y=HDS , color=condition)) +
   geom_boxplot() + theme_bw()  +
  scale_color_colorblind() + theme( text = element_text(size=18)) +
  labs(x = "percent of studies randomized [%]" , y = "centered HDS",
       title = "Randomization effect on HDS\n tested on snRNA-seq Liver Cell Atlas")

ggdens <- ggplot( toPlot , aes( x=random.percent , y=HDS , color=condition)) +
  geom_violin(trim = FALSE) +  #theme_bw()  + 
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1),
               geom = "pointrange",
               position = position_dodge(0.95)
              ) + 
  scale_color_colorblind() + theme( text = element_text(size=18)) +
  labs(x = "percent of studies randomized [%]" , y = "centered HDS",
       title = "Randomization effect on HDS\n tested on snRNA-seq Liver Cell Atlas")

ggsave(filename = 'RandomizationPlotHDS24WeeksNAFLDCohort.png', 
       units = 'px',
       path = '~/Desktop/HepatocyteDamageScore/Results/RandomizationTest/',
       width = 2600, height = 2400)

