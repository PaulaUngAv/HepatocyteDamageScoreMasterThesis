library(ggplot2)
library(ggthemes)
library(Seurat)
library(ggpubr)
library(patchwork)
library("scales")

# start from seurat object that contains only hepatocytes from both conditions
# merged 
# this seurat object was created using a different script called: 
# "SubsetCreationToGenerateListGenesExpressedInHepatocytes.R"
# that was run on cluster server because to heavy for local computer
# 

outputPath <- '~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/'

# if needed
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Quality control 

merged.mouse.lca.hep.nucSeq <- readRDS('~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/countTable_mouseMerged_Hepatocytes_nucSeq.rds')

merged.mouse.lca.hep.nucSeq[["percent.mt"]] <- 
  PercentageFeatureSet(merged.mouse.lca.hep.nucSeq,
                       pattern = "^mt-")

# check dimension before and after next step
dim(merged.mouse.lca.hep.nucSeq)
# genes cells
# 26269 49521
table(merged.mouse.lca.hep.nucSeq@meta.data$orig.ident)
# Nafld   stst 
# 23613  25908 

# keep only sn and sc, we will use sc only as a proxy for the FoldChange() 
# it only if we can split the cells into groups, but we are not really interested 
# in the log fold change between sn and sc but in the other output of the cell
# which is the detection rate of each gene per group, in this case: 
# we get the detection rate of each gene in sn data and that is what need to 
# select which genes we consider as expressed in hepatocytes

merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, 
                                      subset = typeSample == 'nucSeq' )
dim(merged.mouse.lca.hep.nucSeq)
# genes nuclei
# 26269 31604
table(merged.mouse.lca.hep.nucSeq@meta.data$orig.ident)
# nafld  stst
# 20775  10829

merged.mouse.lca.hep.nucSeq <- 
  RenameIdents(object = merged.mouse.lca.hep.nucSeq, 
               'liver_cell_atlas_Nafld' = 'NAFLD')

merged.mouse.lca.hep.nucSeq <- 
  RenameIdents(object = merged.mouse.lca.hep.nucSeq, 
               'liver_cell_atlas_stst' = 'Std. Diet')
#######
# Panel with QC plots
#######
ptemp1 <- VlnPlot(merged.mouse.lca.hep.nucSeq, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt" ), 
        ncol = 3, 
        pt.size = F) & 
  scale_fill_colorblind() & 
  theme( text = element_text(size=16), 
         plot.title = element_text(size = 16),
         axis.text = element_text(size = 10))

ptemp2 <- FeatureScatter(merged.mouse.lca.hep.nucSeq, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt",
                        jitter = TRUE) +  
  scale_color_colorblind() + 
  theme( text = element_text(size=16), legend.position = 'none',
         axis.text = element_text(size = 10))

ptemp2 <- FeatureScatter(subset(merged.mouse.lca.hep.nucSeq, 
                                subset= orig.ident == "liver_cell_atlas_stst"),
                         feature1 = "nCount_RNA", 
                         feature2 = "percent.mt", 
                         group.by="ident",
                         cols = cbbPalette[1],
                         jitter = TRUE) +
  theme( text = element_text(size=16), legend.position = 'none', 
         axis.text = element_text(size = 10)) 

ptemp3 <- FeatureScatter(subset(merged.mouse.lca.hep.nucSeq, 
                      subset= orig.ident == "liver_cell_atlas_Nafld"), 
               feature1 = "nCount_RNA", 
               feature2 = "percent.mt", 
               group.by="orig.ident",
               cols = cbbPalette[2],
               jitter = TRUE) +
  theme( text = element_text(size=16), legend.position = 'none',
         axis.text = element_text(size = 10)) 
  
ptemp4 <- FeatureScatter(merged.mouse.lca.hep.nucSeq, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA", 
                        jitter = TRUE) +  
  scale_color_colorblind() + 
  theme( text = element_text(size=16), legend.position = 'right',
         axis.text = element_text(size = 10)) 

  

qcPlot <- ptemp1 / (ptemp2 |ptemp3 | ptemp4)

ggsave('nucSeqMergedQCPriorToFiltering.png',
       path = outputPath) 

#####
#Filtering after QC
#####

# will use 0.5 mt percentage to be quite stringent but not brutal

# testing everything on dummy
merged.mouse.lca.hep.nucSeq.filtered <- subset(merged.mouse.lca.hep.nucSeq,
                                            nFeature_RNA > 200 &
                                         nFeature_RNA <= 6000 &
                                         percent.mt <= 0.5 & 
                                         nCount_RNA < 30000 )

dim(merged.mouse.lca.hep.nucSeq.filtered)
# genes nuclei
# 26269 20731
table(merged.mouse.lca.hep.nucSeq.filtered@meta.data$orig.ident)
# nafld  stst
# 10485  10246

#####
# Calculate metrics for generating list of genes expressed in hepatocyte
#####
# no need to normalize because I only want to know if gene is detected o not

mergedFilteredMatrix <- GetAssayData(merged.mouse.lca.hep.nucSeq.filtered, 
                                     slot = 'counts')


nRow <- nrow(mergedFilteredMatrix)
#  WARNING: the following command takes 30 minutes to run
res <- sapply(1:nRow, FUN=function(i){
  x <- mergedFilteredMatrix[i,]
  if(i%%100 == 0){print(i)}
  return(sum(x > 0))
})

nCol <- ncol(mergedFilteredMatrix)
metricsHepGenes <- data.frame('Genes' = rownames(mergedFilteredMatrix),
                              'nNuclei' = res, 
                              'percentNuclei' = sapply(1:length(res), 
                                                       FUN = function(i){
                                                         return((res[i]/nCol)*100)
                                                         })
                              )

summary(metricsHepGenes$percentNuclei)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00000  0.02412  0.42931  9.18627  9.92234 99.88906 

gghistogram(data = metricsHepGenes, x = 'percentNuclei', bins= 200) 

finalList <- metricsHepGenes[metricsHepGenes$percentNuclei >= 25, ]
finalList <- finalList[order(finalList$percentNuclei, decreasing = TRUE),]
## save Hepa Genes my method
write.csv(finalList, 
          file = paste0(outputPath,'/GeneListExpressedin25orMoreNuclei.csv' ))
gghistogram(data = finalList, x = 'percentNuclei', bins= 100) + 
  ggtitle('Distribution of nr. genes expressed in 25% or more nuclei')

ggsave('GenesExpressendIn25PercentOrMoreNucleiHistogram.png',
       path = outputPath) 


#### TIM's method 

### get list of genes detected in all experiments

threshold <- ncol(mergedFilteredMatrix)*0.005
expr <- mergedFilteredMatrix[ rowSums( mergedFilteredMatrix > 0 ) > threshold , ]
percentNucleiExpressingGene <- rowSums(mergedFilteredMatrix > 0 )/ ncol(mergedFilteredMatrix)


HepGenes <- data.frame('Genes' = 
                         rownames(expr),
                       'percentNuclei' = 
                         unname(percentNucleiExpressingGene[rownames(expr)]))


write.csv(HepGenes, 
          file = paste0(outputPath,
                   'GeneListTimsMethod.csv' ))

gghistogram(data = HepGenes, x = 'percentNuclei', bins= 100) + 
  ggtitle('Distribution of nr. genes expressed in 0.5% or more nuclei')

ggsave('GenesExpressendIn05PercentOrMoreNucleiHistogramTimsMethod.png',
       path = outputPath) 
