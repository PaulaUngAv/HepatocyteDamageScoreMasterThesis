### Applying HDS on Spatial Transcriptomics ###
library(Seurat)
library(hdf5r)
options(Seurat.object.assay.version = "v5")
library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)
library(gridExtra)
library(AUCell)
library(GSEABase)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(viridis)

# functions
source(file='~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R')
source(file='~/Repositories/cell-damage-score/AUCell_script.r')

options(future.globals.maxSize = 1e+09)
#############################
# 1. Liver Cell Atlas 
# Visium 36 Weeks Mouse Liver
#############################

{
  ###########
  # Load data 
  ###########
  setwd('~/Desktop/HepatocyteDamageScore/SpatialTranscriptomics/LiverCellAtlasVisium/')
  dataPath <- 'GSE192742_36WeeksMouseLiver/'
  plotsPath <- 'GSE192742_36WeeksMouseLiver/Plots/'
  
  # GEO data info: 
  #vstrain: C57Bl/6
  # platform: 10x Visium
  # condition: SD, number of added abs: 0, number of spots: 1279
  # Raw sequencing data for each sample was converted 
  # to matrices of expression counts using the Space Ranger software provided by 10X Genomics (version 1.0).
  # Genome_build: mm10 or hg19
  
  # Sample Info
  # StMouse001Sample001 -> shortfilename: JBO6
  # WDMouse001Sample001 -> shortfilename: JBO9
  # WDMouse002Sample001 -> shortfilename: JBO10
  # WDMouse002Sample002 -> shortfilename: JBO12
  
  # there is another batch of 4 standard diet mouse livers -> GSE192741
  # this are al
  

  annotationNAFLD <- read.csv(file = paste0(dataPath,"annot_mouseNafldVisium.csv"))

 
  spatialSamples <- list(
    Load10X_Spatial(data.dir = paste0(dataPath,'StSt/Mouse001Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5",
                    assay = "LiverCellAtlasSpatial",
                    slice = "StMouse001Sample001",
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "StSt/Mouse001Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse001Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5", 
                    assay = "LiverCellAtlasSpatial", 
                    slice = "WDMouse001Sample001", 
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse001Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse002Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5",
                    assay = "LiverCellAtlasSpatial",
                    slice = "WDMouse002Sample001",
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse002Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse002Sample002/'),
                    filename = "filtered_feature_bc_matrix.h5", 
                    assay = "LiverCellAtlasSpatial", 
                    slice = "WDMouse002Sample002", 
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse002Sample002/")))
    )
  
  sampleNames <- c("StSpatialSeuratM1S1",
                   "WDSpatialSeuratM1S1",
                   "WDSpatialSeuratM2S1", 
                   "WDSpatialSeuratM2S2")
  
  names(spatialSamples) <- sampleNames

  sampleCoordinates <- list(
    spatialSamples[[1]]@images$StMouse001Sample001@coordinates,
    spatialSamples[[2]]@images$WDMouse001Sample001@coordinates,
    spatialSamples[[3]]@images$WDMouse002Sample001@coordinates,
    spatialSamples[[4]]@images$WDMouse002Sample002@coordinates
    )
  
  ### Functions from Liver Cell Atlas Git Hub
  ##### Function drawVlnPlot
  
  drawVlnPlot <- function(toPlot, colsToColor){
    toPlot$staticNr=1
    
    toPlot<-toPlot[order(toPlot[,colsToColor[1]]),]
    p_nGene <- ggplot(toPlot, aes(staticNr, nFeature_LiverCellAtlasSpatial)) + 
      geom_violin(fill="gray80") + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5, size=0.8) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic() +
      theme(legend.position = "none", axis.title.x = element_blank())
    
    toPlot<-toPlot[order(toPlot[,colsToColor[2]]),]
    p_nUMI <- ggplot(toPlot, aes(staticNr, nCount_LiverCellAtlasSpatial)) + 
      geom_violin(fill="gray80") + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5, size=0.8) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic() +
      theme(legend.position = "none", axis.title.x = element_blank())
    
    toPlot<-toPlot[order(toPlot[,colsToColor[3]]),]
    p_mito <- ggplot(toPlot, aes(staticNr, percent.mt)) + 
      geom_violin(fill="gray80") + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5, size=0.8) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic() +
      theme(legend.position = "none", axis.title.x = element_blank())
    
    return(grid.arrange(p_nGene,p_nUMI,p_mito,ncol=3))
  }
  
  
  ##### Function drawVlnPlot_split
  drawVlnPlot_split <- function(toPlot, colsToColor){
    toPlot$staticNr=1
    
    columnNr<-which(colnames(toPlot)==colsToColor[1])
    toPlot<-toPlot[order(toPlot[,columnNr]),]
    p_nGene <- ggplot(toPlot, aes(orig.ident, nFeature_LiverCellAtlasSpatial, fill=orig.ident)) + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5, size=0.8) +
      geom_violin() + 
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      theme_classic() +
      theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
            axis.text.y = element_text(size = 8), axis.title.x = element_blank())
    
    columnNr<-which(colnames(toPlot)==colsToColor[2])
    toPlot<-toPlot[order(toPlot[,columnNr]),]
    p_nUMI <- ggplot(toPlot, aes(orig.ident, nCount_LiverCellAtlasSpatial, fill=orig.ident)) + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5, size=0.8) +
      geom_violin() + 
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      theme_classic() +
      theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
            axis.text.y = element_text(size = 8), axis.title.x = element_blank())
    
    columnNr<-which(colnames(toPlot)==colsToColor[3])
    toPlot<-toPlot[order(toPlot[,columnNr]),]
    p_mito <- ggplot(toPlot, aes(orig.ident, percent.mt, fill=orig.ident)) + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5, size=0.8) +
      geom_violin() +
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      theme_classic() +
      theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
            axis.text.y = element_text(size = 8), axis.title.x = element_blank())
    
    thePlot<-grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
    return(thePlot)
  }
  
  ########################################  
  #### Check for outlier sports (from liver cell atlas git hub)
  ########################################
  
  ##### Sample: StSpatialSeurateM1S1 #####
  
  # p1 <- SpatialDimPlot(spatialSamples[[1]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[1]])+
  #   geom_point(size=0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  # 
  # grid.arrange(p1, p2, ncol=2)
  # 
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #        file='/PATH/TO/results/image_sample1.png', 
  #        width=10, height = 5)
  # 
  # Select outlier spots
  badSpots <- sampleCoordinates[[1]] %>% mutate('spot' = rownames(.)) %>% 
    dplyr::filter(imagecol < 5000 | imagerow < 4500) %>%
    pull(spot)
  
  # Visualize bad spots
  sampleCoordinates[[1]]$badSpot <- FALSE
  sampleCoordinates[[1]][badSpots,'badSpot'] <- TRUE
  # p2 <- ggplot(aes(x = imagecol, y = imagerow, color = badSpot), 
  #              data = sampleCoordinates[[1]]) + 
  #   geom_point(size=0.5) +
  #   theme_classic() + 
  #   scale_y_reverse()
  # 
  # grid.arrange( p1, p2, ncol = 2)
  
  # Remove bad spots
  spatialSamples[[1]] <- spatialSamples[[1]][,setdiff(colnames(spatialSamples[[1]]),badSpots)]
  dim(spatialSamples[[1]])
  # 31053  1274
  
  # p3 <- SpatialDimPlot(spatialSamples[[1]])
  # ggsave(grid.arrange(p1,p2,p3,ncol=3), file='/PATH/TO/results/image_sample2.png', 
  #       width=15, height = 5)
  
  
  ##### Sample: WDSpatialSeuratM1S1 #####

  # p1 <- SpatialDimPlot(spatialSamples[[2]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[2]])+
  #   geom_point(size = 0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  
  # grid.arrange(p1, p2, ncol=2)
  
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #       file='/PATH/TO/results/image_sample1.png', 
  #       width=10, height = 5)
  
  # No outlier spots
  
  #### Sample: WDSpatialSeuratM2S1 ####
  # 
  # p1 <- SpatialDimPlot(spatialSamples[[3]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[3]])+
  #   geom_point(size = 0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  # 
  # grid.arrange(p1, p2, ncol = 2)
  
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #       file='/PATH/TO/results/image_sample1.png', 
  #       width=10, height = 5)
  
  # Select outlier spots
  badSpots <- sampleCoordinates[[3]] %>% mutate('spot' = rownames(.)) %>% 
    dplyr::filter(imagerow > 16300) %>%
    pull(spot)
  
  # Visualize bad spots
  sampleCoordinates[[3]]$badSpot <- FALSE
  sampleCoordinates[[3]][badSpots,'badSpot'] <- TRUE
  # p2 <- ggplot(aes(x = imagecol, y = imagerow, color = badSpot), 
  #              data = sampleCoordinates[[3]]) + 
  #   geom_point(size=0.5) +
  #   theme_classic() + 
  #   scale_y_reverse()
  # 
  # grid.arrange( p1, p2, ncol = 2)
  
  # Remove bad spots
  spatialSamples[[3]] <- spatialSamples[[3]][,setdiff(colnames(spatialSamples[[3]]),badSpots)]
  dim(spatialSamples[[3]])
  
  # 31053  2074
  
  #p3 <- SpatialDimPlot(spatialSamples[[3]])
  # ggsave(grid.arrange(p1,p2,p3,ncol=3), file='/PATH/TO/results/image_sample2.png', 
  #       width=15, height = 5)
  

  
  #### Sample: WDSpatialSeuratM2S2 ####
  
  # p1 <- SpatialDimPlot(spatialSamples[[4]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[4]])+
  #   geom_point(size = 0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  # 
  # grid.arrange(p1, p2, ncol = 2)
  # 
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #       file='/PATH/TO/results/image_sample1.png', 
  #       width=10, height = 5)
  
  # Select outlier spots
  badSpots <- sampleCoordinates[[4]] %>% mutate('spot' = rownames(.)) %>% 
    dplyr::filter(imagecol < 5000 | imagerow > 16000) %>%
    pull(spot)
  
  # Visualize bad spots
  sampleCoordinates[[4]]$badSpot <- FALSE
  sampleCoordinates[[4]][badSpots,'badSpot'] <- TRUE
  # p2 <- ggplot(aes(x = imagecol, y = imagerow, color = badSpot), 
  #              data = sampleCoordinates[[4]]) + 
  #   geom_point(size = 0.5) +
  #   theme_classic() + 
  #   scale_y_reverse()
  # 
  # grid.arrange( p1, p2, ncol = 2)
  # 
   # Remove bad spots
  spatialSamples[[4]] <- 
    spatialSamples[[4]][, setdiff(colnames(spatialSamples[[4]]),badSpots)]
  
  dim(spatialSamples[[4]])
  # 31053  2252
  
  #p3 <- SpatialDimPlot(spatialSamples[[4]])
  # ggsave(grid.arrange(p1,p2,p3,ncol=3), file='/PATH/TO/results/image_sample2.png', 
  #       width=15, height = 5)
  
  
  ########################################
  # Merge 
  ########################################
  
  spatialSamples[[1]]@meta.data$orig.ident <- sampleNames[1]
  spatialSamples[[2]]@meta.data$orig.ident <- sampleNames[2]
  spatialSamples[[3]]@meta.data$orig.ident <- sampleNames[3]
  spatialSamples[[4]]@meta.data$orig.ident <- sampleNames[4]
  


  spatialObjectMerged <- merge(spatialSamples[[1]], y = c(spatialSamples[[2]], 
                                                      spatialSamples[[3]],
                                                      spatialSamples[[4]]), 
                               add.cell.ids = sampleNames, merge.data = TRUE)
  
  # checking merge
  dim(spatialObjectMerged)
  head(spatialObjectMerged@meta.data)
  tail(spatialObjectMerged@meta.data)
  table(spatialObjectMerged@meta.data$orig.ident)
  
  # StSpatialSeuratM1S1 WDSpatialSeuratM1S1 WDSpatialSeuratM2S1 WDSpatialSeuratM2S2 
  # 1274                2348                2073                2252 
  
  
  #########################################
  # Filtering
  #########################################
  
  spatialObjectMerged[["percent.mt"]] <- 
    PercentageFeatureSet(spatialObjectMerged, pattern = "^mt-")
  spatialObjectMerged
  
  # Select spots to be removed: kept same thresholds as Liver Cell Atlas git hub
  
  toPlot <- spatialObjectMerged@meta.data
  toPlot$badCell1 <- FALSE
  toPlot[ toPlot$nFeature_LiverCellAtlasSpatial < 1000, 'badCell1' ] <- TRUE
  toPlot$badCell2 <- FALSE
  toPlot[toPlot$nCount_LiverCellAtlasSpatial < 5000, 'badCell2'] <- TRUE
  toPlot[toPlot$nCount_LiverCellAtlasSpatial> 60000, 'badCell2'] <- TRUE
  toPlot$badCell3 <- FALSE
  toPlot[toPlot$percent.mt > 15, 'badCell3'] <- TRUE
  
  # Here I use functions from Liver Cell Atlas GitHub
  
  p1 <- drawVlnPlot(toPlot,colsToColor = c('badCell1','badCell2','badCell3'))
  print(p1)
  ggsave(p1,file='SpotsOverThresholdSCT.png')
  

  # Do filtering
  spatialObjectMergedFiltered <-subset(spatialObjectMerged, 
                         nFeature_LiverCellAtlasSpatial > 1000 & 
                           nCount_LiverCellAtlasSpatial > 5000 & 
                           nCount_LiverCellAtlasSpatial < 60000 & 
                           percent.mt < 15)
  
  dim(spatialObjectMergedFiltered)
  # From   31053  7947 to 31053  7828
  
  ########################################
  # Normalization with SCT 
  ########################################

  ## the reason for the normalization pre filtering is that if we normalize
  ## post merging we normalize all sections together 
  
  spatialObjectMergedFiltered <- SCTransform(spatialObjectMergedFiltered , 
                                             verbose = TRUE, 
                                             assay = "LiverCellAtlasSpatial")
  remove(spatialObjectMerged)
  
  # this would be the alternative: to filter and transform each sample individually
  # I will compare the UMAPs of both methods to check if there is an effect of 
  # normalizing the data when merged or not-merged
  
  # do filtering and SCT normalization with separate seurat objects in case
  # we need them 
  
  i = 1
  while(i<=4){
    spatialSamples[[i]][["percent.mt"]] <- 
      PercentageFeatureSet(spatialSamples[[i]], pattern = "^mt-")
    spatialSamples[[i]] <- subset(spatialSamples[[i]], 
                                  nFeature_LiverCellAtlasSpatial > 1000 & 
                                    nCount_LiverCellAtlasSpatial > 5000 & 
                                    nCount_LiverCellAtlasSpatial < 60000 & 
                                    percent.mt < 15)
    spatialSamples[[i]] <- SCTransform(spatialSamples[[i]], 
                                       verbose = TRUE, 
                                       assay = "LiverCellAtlasSpatial")
    
    i <- i + 1
    
  }

  
  ###########################
  # Add annotation info to meta data (from GEO) of merged and filtered object
  ###########################
  
  # the annotation file from GEO includes a zonation parameter, a zonation group 
  # and cluster identity for each cell 
  # here I noticed that there are less spots annotated than there are spots in 
  # spots in the data set, so probably they got removed
  
  # reformat meta data
  spatialObjectMergedFiltered@meta.data$barcode <- 
    sub( ".*_", "", rownames(spatialObjectMergedFiltered@meta.data))
  
  # barcodes are repeated between samples! i.e. not unique, so need sample
  # names added to barcode to differentiate between spots
  
  spatialObjectMergedFiltered@meta.data$sample_barcode <- 
    rownames(spatialObjectMergedFiltered@meta.data)
  
  #do this for each samples
  #length(which(annotationNAFLD$sample == 'JBO6'))
  # [1] 1121
  #length(which(spatialObjectMergedFiltered$orig.ident == 'StSpatialSeuratM1S1'))
  # [1] 1244
  annotationNAFLD$sample_barcode <- NA
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO6'),'sample_barcode'] <-
    paste0('StSpatialSeuratM1S1_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO6'),'spot'])
  
  #length(which(annotationNAFLD$sample == 'JBO9'))
  # 2002
  # length(which(spatialObjectMergedFiltered$orig.ident == 'WDSpatialSeuratM1S1'))
  # 2340
  
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO9'),'sample_barcode'] <-
    paste0('WDSpatialSeuratM1S1_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO9'),'spot'])
  
  #length(which(annotationNAFLD$sample == 'JBO10'))
  # 1780
  #length(which(spatialObjectMergedFiltered$orig.ident == 'WDSpatialSeuratM2S1'))
  # 2020
  
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO10'),'sample_barcode'] <-
    paste0('WDSpatialSeuratM2S1_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO10'),'spot'])
  
  #length(which(annotationNAFLD$sample == 'JBO12'))
  # 1768
  #length(which(spatialObjectMergedFiltered$orig.ident == 'WDSpatialSeuratM2S2'))
  # 2224
  
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO12'),'sample_barcode'] <-
    paste0('WDSpatialSeuratM2S2_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO12'),'spot'])
  
  annotationNAFLD$sample_barcode <- gsub(pattern = '*_1|_2|_3|_4', 
                                         x = annotationNAFLD$sample_barcode, 
                                         replacement = '')
  
  # match indices to add info to meta data
  spatialObjectMergedFiltered@meta.data$annot_cluster <- NA
  spatialObjectMergedFiltered@meta.data$zonation <- NA
  spatialObjectMergedFiltered@meta.data$zonationGroup <- NA
  
  copyMetadata <- spatialObjectMergedFiltered@meta.data
  i = 1
  while(i <= length(annotationNAFLD$sample_barcode)){
    if(sum(copyMetadata$sample_barcode == annotationNAFLD$sample_barcode[i]) == 1){
      copyMetadata$annot_cluster[
        copyMetadata$sample_barcode == 
          annotationNAFLD$sample_barcode[i]] <- annotationNAFLD$cluster[i]
      copyMetadata$zonation[
        copyMetadata$sample_barcode == 
          annotationNAFLD$sample_barcode[i]] <- annotationNAFLD$zonation[i]
      copyMetadata$zonationGroup[
        copyMetadata$sample_barcode == 
          annotationNAFLD$sample_barcode[i]] <- annotationNAFLD$zonationGroup[i]
    }else{
     print('not there!')
    print(i)
    }
    i = i +1
  }
  
  # add column to meta data 
  if(identical(spatialObjectMergedFiltered@meta.data$sample_barcode, 
  copyMetadata$sample_barcode)){
    spatialObjectMergedFiltered@meta.data$annot_cluster <- 
      as.factor(copyMetadata$annot_cluster)
    spatialObjectMergedFiltered@meta.data$zonation <- 
      copyMetadata$zonation
    spatialObjectMergedFiltered@meta.data$zonationGroup <- 
      as.factor(copyMetadata$zonationGroup)
  }
  


  ###########################
  # UMAP
  ###########################
  
  # UMAP of merged and then transformed "spatialObjectMergedFiltered"
  
  # alternative 1: merged and then normalized
  
  DefaultAssay(spatialObjectMergedFiltered) <- "SCT"
  
  # this step is done in the github
  ### !!! not sure about this step 
  VariableFeatures(spatialObjectMergedFiltered ) <- c(VariableFeatures(spatialSamples[[1]]), 
                                        VariableFeatures(spatialSamples[[2]]),
                                        VariableFeatures(spatialSamples[[3]]),
                                        VariableFeatures(spatialSamples[[4]]))
  
  spatialObjectMergedFiltered  <- RunPCA(spatialObjectMergedFiltered, 
                                         verbose = FALSE)
  
  spatialObjectMergedFiltered  <- FindNeighbors(spatialObjectMergedFiltered , 
                                                dims = 1:30)
  
  spatialObjectMergedFiltered  <- FindClusters(spatialObjectMergedFiltered , 
                                               resolution = 0.6)
  
  spatialObjectMergedFiltered  <- RunUMAP(spatialObjectMergedFiltered , dims = 1:30)
  
  
  # calculates percetage based on total mt in all samples merged (right?)
  spatialObjectMergedFiltered [["percent.mt_merged"]] <- 
    PercentageFeatureSet(spatialObjectMergedFiltered, pattern = "^mt-")
  
  ## RBC stands for 'red blood cells', hemoglobin (different ones) expressing cells
  # I would have to look at this in paper 
  spatialObjectMergedFiltered [["percent.rbc_merged"]] <- 
    PercentageFeatureSet(spatialObjectMergedFiltered, pattern = "Hb[ab]-")
  
  plotMito <- FeaturePlot(spatialObjectMergedFiltered, c("percent.mt"))
  plotMitoMerged <- FeaturePlot(spatialObjectMergedFiltered, c("percent.mt_merged"))
  # plots look similar 
  
  plotRBC <- FeaturePlot(spatialObjectMergedFiltered, c("percent.rbc_merged"))
  
  ########################################
  # Plot on UMAP (all of this just like in github from Liver Cell Atlas)
  ########################################
  
  ### Plot clusters
  umapPlot <- DimPlot(spatialObjectMergedFiltered, 
                      reduction = "umap",
                      order = c(4,6),
                      label = T)
  ggsave(filename = 'UMAPFilteredMergedSCT.png',
        path = plotsPath)
  
  umapPlotSplit <- DimPlot(spatialObjectMergedFiltered, 
                           reduction = "umap",
                           order = c(4,6),
                           group.by = 'orig.ident')
  
  pcaPlot <- DimPlot(spatialObjectMergedFiltered, 
          reduction = "pca",
          pt.size = 0.1, 
          group.by = 'orig.ident',
          shuffle = TRUE) + scale_color_colorblind() +
    theme(legend.justification=c(-0.1,1), legend.position=c(0,0.95)) + 
    labs(title = NULL)
    
  
  ggsave(filename = 'PCAFilteredMergedSCT.png',
         path = plotsPath)
  
  
  # plotting clusters predominant in WD samples last to see if they are also 
  # present in left cluster -> it seems like there are not clsuters 4 or 5 in left clsuter
  
  DimPlot(spatialObjectMergedFiltered, 
          reduction = "umap",
          order = c(4,6),
          pt.size = 0.7
          )
  
  PCAPlot(spatialObjectMergedFiltered,dims = c(3,4))
  
  grid.arrange(umapPlot, 
               umapPlotSplit, 
               nrow=2)
  
  ggsave(filename = 'UMAPandSampleIdentitiesFilteredMergedSCT.png',
         path = plotsPath)
  
  ### Plot mitochondrial genes percentage on UMAP
  # when explaining plots, dont forget to mention cut-offs
  plotMito <- FeaturePlot(spatialObjectMergedFiltered,'percent.mt_merged',
                        min.cutoff = 'q2',  
                        max.cutoff = 'q98')
  
  grid.arrange(umapPlot, umapPlotSplit, plotMito, ncol=3)
  
  ### Plot RBC
  # when explaining plots, dont forget to mention cut-offs
  plotRBC <- FeaturePlot(spatialObjectMergedFiltered, 
                         'percent.rbc_merged', 
                         min.cutoff = 'q2', 
                         max.cutoff = 'q98')
  
  grid.arrange(umapPlot, umapPlotSplit, plotMito, plotRBC, ncol=4)
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, plotMito, 
                      plotRBC, ncol=2, nrow=2),
         file=paste0(dataPath, '/spatialMergedFilteredUMAPall.png'), 
         width = 24, height = 6)
  
  
  ### Plot some genes
  FeaturePlot(spatialObjectMergedFiltered, c("Glul", "Cyp2f2"), 
              min.cutoff = 'q2', max.cutoff = 'q98')
  
  
  ### Plot clusters
  SpatialDimPlot(spatialObjectMergedFiltered,
                 idents)

  
  #grid.arrange(SpatialDimPlot(spatialObjectMergedFiltered), ncol = 2,
  #             nrow = 2)
  SpatialDimPlot(spatialObjectMergedFiltered,label = TRUE, label.size = 3)
  
  
  #Plot certain clusters (this are numbers on git hub, probably other numbers in here)
  # clustersOI<-c(6,9,5,8,3,11)
  # 
  # SpatialDimPlot(spatialObjectMergedFiltered, 
  #                cells.highlight = CellsByIdentities(
  #                  object = spatialObjectMergedFiltered, 
  #                  idents = clustersOI), 
  #                facet.highlight = TRUE, 
  #                images = 'StMouse001Sample001', 
  #                ncol=6)
  # 
  #  
  # SpatialDimPlot(spatialObjectMergedFiltered, 
  #                cells.highlight = CellsByIdentities(
  #                  object = spatialObjectMergedFiltered, 
  #                  idents = clustersOI), 
  #                facet.highlight = TRUE, 
  #                images = 'WDMouse001Sample001',
  #                ncol=6)
  # 
  # SpatialDimPlot(spatialObjectMergedFiltered, 
  #                cells.highlight = CellsByIdentities(
  #                  object = spatialObjectMergedFiltered, 
  #                  idents = clustersOI), 
  #                facet.highlight = TRUE, 
  #                images = 'WDMouse002Sample001', 
  #                ncol=6)
  # 
  # ### Plot mito
  # SpatialFeaturePlot(spatialObjectMergedFiltered, 
  #                    features = "percent.mt") 
  
  
  ### Plot genes (why theses genes?)
  SpatialFeaturePlot(spatialObjectMergedFiltered, 
                     features = c("Glul", "Cyp2f2"))
  
  ####################
  # Apply HDS
  ####################
  
  # Genes expressed in hepatocytes
  geneFilter <- read.csv(
    '~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')
  geneFilter <- geneFilter$Genes
  
  # Load HUDS (Hepatocyte Universal Damage Score)
  geneListHUDS <- read.csv(
    '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')
  
  # extract expression matrix 
  exprMatrix <- GetAssayData(spatialObjectMergedFiltered, slot = "counts",
                             assay = "SCT")
  
  HDSmergedfilteredSCT <- DS_calc.func( exprMatrices = exprMatrix, 
                           ceilThrsh = 0.05 , 
                           DSignature = geneListHUDS, 
                           ntop =42 , 
                           useOrder = 'mean_rank')
  
  HDSmergedfilteredSCT <- data.frame('cellID' = names(HDSmergedfilteredSCT),
                                     'HDS_SCT' = HDSmergedfilteredSCT )
  
  # add HDS to SueratObjectMetadata for plotting 
  identical(row.names(spatialObjectMergedFiltered@meta.data), 
            HDSmergedfilteredSCT$cellID)
  # TRUE
  
  spatialObjectMergedFiltered@meta.data$HDS_SCT <- HDSmergedfilteredSCT$HDS_SCT
  
  ########################################
  # Plot HDS
  ########################################

  # format and center data for box-plots

  HDSmergedfilteredSCT$orig.ident <- 
    spatialObjectMergedFiltered@meta.data$orig.ident
  HDSmergedfilteredSCT$centeredHDS <- scale(HDSmergedfilteredSCT$HDS_SCT, 
                                            scale = FALSE)
  HDSmergedfilteredSCT$zonationGroup <- 
    spatialObjectMergedFiltered@meta.data$zonationGroup
  
  HDSmergedfilteredSCT$zonation <- 
    spatialObjectMergedFiltered@meta.data$zonation
  
  HDSmergedfilteredSCT$diet <- NA
  HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == 'StSpatialSeuratM1S1', 'diet'] <- 
    'Standard Diet' 
  HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident != 'StSpatialSeuratM1S1', 'diet'] <- 
    'Western Diet'
  HDSmergedfilteredSCT$diet <- factor(HDSmergedfilteredSCT$diet,
                                      ordered = TRUE, 
                                      levels = c("Standard Diet", "Western Diet"))

  
  # Violin Plot of HDS distribution
  
  CenteredHDSViolinPlot <- ggplot(data = HDSmergedfilteredSCT,
                           aes(x = orig.ident, y = centeredHDS)) +
    geom_violin() + xlab("samples") + ylab("centered HDS") +
    scale_x_discrete(labels = c("Standard Diet M1S1", "Western Diet M1S1",
                                "Western Diet M2S1", "Western Diet M2S2")) +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12))  + 
    geom_boxplot(width=0.1, outlier.size = 0) 
    
    
  ggsave(CenteredHDSViolinPlot, 
         path = plotsPath,
         file= 'CenteredHDSViolinPlotMergedFilteredSCT.png')
  
  # let's try not centered

  HDSViolinPlot <- ggplot(data = HDSmergedfilteredSCT,
                          aes(x = orig.ident, y = HDS_SCT,
                              color = diet)) +
    geom_violin() + xlab("samples") + ylab("HDS") +
    scale_x_discrete(labels = c("St.D. Mouse 1 S1", "W.D. Mouse 1 S1",
                                "W.D. Mouse 2 S1", "W.D. Mouse 2 S2")) +
    theme(axis.text.x = element_text(size=14, angle = 45),
          axis.text.y = element_text(size=14),
          legend.position = 'none')  + 
    scale_color_colorblind() +
    geom_boxplot(width = 0.1, outlier.size = 0) 
  
  ggsave(HDSViolinPlot, 
         path = plotsPath,
         file= 'HDSViolinPlotMergedFilteredSCT.png')
 
 # DONT THINKS I'LL DO THIS: Plot of spatial distribution of HDS centered around zero
  
  spatialObjectMergedFiltered@meta.data$centeredHDS_SCT <- 
    HDSmergedfilteredSCT$centeredHDS
  
  # add diet column to meta deta
  
  if(identical(HDSmergedfilteredSCT$cellID, 
               rownames(spatialObjectMergedFiltered@meta.data))){
    spatialObjectMergedFiltered@meta.data$diet <- HDSmergedfilteredSCT$diet
  }
  
  # Plot HDS on spots 
  
  tempPlot <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                                 features = "HDS_SCT"
  ) &
    scale_fill_gradientn(
      colours =  c("#FFFFCC","#a1dab4","#41b6c4","#2c7fb8","#253494"),
      values = NULL,
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
     limits = c(-0.47387,0.05601)
    ) 
  # try viridis 
  # I like it a lot!!
  tempPlot <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                                 features = "HDS_SCT"
  ) &
    scale_fill_gradientn(
      colours = turbo(n = 100),
      values = NULL,
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = c(-0.47387,0.05601)
    ) 
  
  # plots for thesis:
  
  
  p1 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS_SCT", 
                           images = 'StMouse001Sample001', 
                           pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  turbo(n = 100),
                          aesthetics = "fill", 
                          limits = c(-0.47387,0.05601),
                          name = "HDS") +
    theme(legend.position = 'none') +
    ggtitle("Standard Diet\nMouse 1 Sample 1")
  
  
  p2 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS_SCT",
                           images = 'WDMouse001Sample001',
                           pt.size.factor = 1.8) + 
    scale_fill_gradientn( colours =  turbo(n = 100),
                          aesthetics = "fill", 
                          limits = c(-0.47387,0.05601),
                          name = "HDS")+
    theme(legend.position = 'none'
          # plot.title = element_text(size = 20)
    ) +
    ggtitle("Western Diet\nMouse 1 Sample 1")
  
  p3 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS_SCT",
                           images = 'WDMouse002Sample001',
                           pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  turbo(n = 100),
                          aesthetics = "fill", 
                          limits = c(-0.47387,0.05601),
                          name = "HDS") +
    theme(legend.position = 'none') +
    ggtitle("Western Diet\nMouse 2 Sample 1")
  
  p4 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS_SCT",
                           images = 'WDMouse002Sample002',
                           pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  turbo(n = 100),
                          aesthetics = "fill", 
                          limits = c(-0.47387,0.05601),
                          name = "HDS") +
    theme(legend.position = 'right'
          #legend.key.size = unit(1.5, 'cm'),
          #legend.title = element_text(size = 20),
          #plot.title = element_text(size = 20)
    ) +
    ggtitle("Western Diet\nMouse 2 Sample 2")
  
  plotJustImages <- wrap_plots(SpatialPlot(spatialObjectMergedFiltered, 
                                           alpha = 0)  & 
                                 theme(legend.position = "none") & ggtitle(''))
  
  plotMito <- wrap_plots(SpatialFeaturePlot(spatialObjectMergedFiltered,
                                            images = 'StMouse001Sample001',
                                            features = "percent.mt", 
                                            pt.size.factor = 2)  + 
                           scale_fill_gradientn( colours = viridis(n = 100),
                                                 aesthetics = "fill", 
                                                 limits = c(0.4342,
                                                            12.1084))  
                         + ggtitle('') + theme(legend.position = 'none') | 
                           
                           SpatialFeaturePlot(spatialObjectMergedFiltered,
                                              images = 'WDMouse001Sample001',
                                              features = "percent.mt", 
                                              pt.size.factor = 1.8
                           )  +
                           scale_fill_gradientn( colours = viridis(n = 100),
                                                 aesthetics = "fill", 
                                                 limits = c(0.4342,
                                                            12.1084))  + 
                           ggtitle('') + theme(legend.position = 'none') |
                           
                           SpatialFeaturePlot(spatialObjectMergedFiltered,
                                              images = 'WDMouse002Sample001',
                                              features = "percent.mt", 
                                              pt.size.factor = 1.8)  + 
                           scale_fill_gradientn( colours = viridis(n = 100),
                                                 aesthetics = "fill", 
                                                 limits = c(0.4342,
                                                            12.1084)) +
                           ggtitle('') + theme(legend.position = 'none') |
                           
                           SpatialFeaturePlot(spatialObjectMergedFiltered,
                                              images = 'WDMouse002Sample002',
                                              features = "percent.mt", 
                                              pt.size.factor = 1.8)  + 
                           scale_fill_gradientn( colours = viridis(n = 100),
                                                 aesthetics = "fill", 
                                                 limits = c(0.4342,
                                                            12.1084), 
                                                 name = "% mit. RNA") + 
                           ggtitle('') + theme(legend.position = 'right')
  ) 
  
  
  # combine HDS on tissue, tissue without spots and box-plot of HDS distribution
  
  layer1plot <- wrap_plots((p1 | p2 | p3 | p4))
  panelplot <-  layer1plot / plotJustImages/ plotMito 
  panelplot
  ggsave('HDS_percMitonTissuePlusTissueSCTFilteredMerged.png',
         path = plotsPath)
  
  layer1plot / plotJustImages
  ggsave('HDSTissueSCTFilteredMerged.png',
         path = plotsPath)
  
  ## plot annotation clusters 
  
  pp1ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'StMouse001Sample001',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom', legend.title = element_text('cluster ID'))
  pp2ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'WDMouse001Sample001',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom',  legend.title = element_text('cluster ID'))
  pp3ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'WDMouse002Sample001',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom', legend.title = element_text('cluster ID'))
  pp4ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'WDMouse002Sample002',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom', legend.title = element_text('cluster ID'))
  
  
  ClusterIdsTissue <- plotMito / 
    wrap_plots((pp1ClusterIdents | pp2ClusterIdents | pp3ClusterIdents | pp4ClusterIdents))
  
  ggsave('ClusterIdsPlusMitoPercNnTissue.png', path = plotsPath)
  
  saveRDS(spatialObjectMergedFiltered,file = '~/Desktop/HepatocyteDamageScore/Results/Spatial/LiverCellAtlas/HDSOnSpatialTranscriptomicsSeuratObject.rds')
  saveRDS(HDSmergedfilteredSCT, file = '~/Desktop/HepatocyteDamageScore/Results/Spatial/LiverCellAtlas/HDSOnSpatialTranscriptomicsDataFrameHDS.rds')
  ### OTHER PLOTS 

  # plot spots that have HDS higher than threshold

  tempPlot <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                                 features = "centeredHDS_SCT"
                                 ) &
    scale_fill_gradientn(
      colours =  c("#FFFFCC","#a1dab4","#41b6c4","#2c7fb8","#253494"),
      values = NULL,
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = c(-0.34,0.32)
    ) 
  
  # test plot all samples
  SpatialFeaturePlot(spatialObjectMergedFiltered,
                     features = "centeredHDS_SCT") &
    scale_fill_viridis(limits = c(-0.34,0.23),
                       direction = -1, 
                       name = 'centered HDS')
  
  # test plot only 
  subsetForPlot <- subset(x =spatialObjectMergedFiltered, 
                          subset = orig.ident != "StSpatialSeuratM1S1")
  
  SpatialFeaturePlot(subsetForPlot,
                     features = "centeredHDS_SCT") +
    scale_fill_viridis(limits = c(-0.1,0.23),
                       direction = -1, 
                       name = 'centered HDS')
  
  # Plots for panel
  
  p1 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "centeredHDS_SCT", 
                           images = 'StMouse001Sample001', 
                            pt.size.factor = 2.3) + 
    scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                       "#2c7fb8","#253494"),
                          aesthetics = "fill", 
                          limits = c(-0.34,0.32),
                          name = "centered HDS") +
    theme(legend.position = 'none'
          # plot.title = element_text(size = 20)
          ) +
    ggtitle("Standard Diet\nMouse 1 Sample 1")
    
  
  p2 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "centeredHDS_SCT",
                           images = 'WDMouse001Sample001',
                           pt.size.factor = 1.8) + 
    scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                       "#2c7fb8","#253494"),
                          aesthetics = "fill", 
                          limits = c(-0.34,0.32),
                          name = "centered HDS") +
    theme(legend.position = 'none'
          # plot.title = element_text(size = 20)
          ) +
    ggtitle("Western Diet\nMouse 1 Sample 1")
  
  p3 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                          features = "centeredHDS_SCT",
                          images = 'WDMouse002Sample001',
                          pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                       "#2c7fb8","#253494"),
                          aesthetics = "fill", 
                          limits = c(-0.34,0.32),
                          name = "centered HDS") +
    theme(legend.position = 'none'
          # plot.title = element_text(size = 20)
          ) +
    ggtitle("Western Diet\nMouse 2 Sample 1")
  
  p4 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                          features = "centeredHDS_SCT",
                          images = 'WDMouse002Sample002',
                          pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                       "#2c7fb8","#253494"),
                          aesthetics = "fill", 
                          limits = c(-0.34,0.32),
                          name = "centered HDS") +
    theme(legend.position = 'right'
          #legend.key.size = unit(1.5, 'cm'),
          #legend.title = element_text(size = 20),
          #plot.title = element_text(size = 20)
          ) +
    ggtitle("Western Diet\nMouse 2 Sample 2")
  
  plotJustImages <- SpatialPlot(spatialObjectMergedFiltered, 
                                alpha = 0)  & 
    theme(legend.position = "none") & ggtitle('')
  
  plotMito <- wrap_plots(SpatialFeaturePlot(spatialObjectMergedFiltered,
                                            images = 'StMouse001Sample001',
                                            features = "percent.mt", 
                                            pt.size.factor = 2)  + 
                           scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                                              "#2c7fb8","#253494"),
                                                 aesthetics = "fill", 
                                                 limits = c(0.468007,
                                                            9.691789)) 
                         + ggtitle('') + theme(legend.position = 'none') | 
                           
                           SpatialFeaturePlot(spatialObjectMergedFiltered,
                                              images = 'WDMouse001Sample001',
                                              features = "percent.mt", 
                                              pt.size.factor = 1.8
                                              )  +
                           scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                                              "#2c7fb8","#253494"),
                                                 aesthetics = "fill", 
                                                 limits = c(0.468007,
                                                            9.691789)) + 
                           ggtitle('') + theme(legend.position = 'none') |
                           
                           SpatialFeaturePlot(spatialObjectMergedFiltered,
                                              images = 'WDMouse002Sample001',
                                              features = "percent.mt", 
                                              pt.size.factor = 1.8)  + 
                           scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                                              "#2c7fb8","#253494"),
                                                 aesthetics = "fill", 
                                                 limits = c(0.468007,
                                                            9.691789)) +
                           ggtitle('') + theme(legend.position = 'none') |
                           
                           SpatialFeaturePlot(spatialObjectMergedFiltered,
                                              images = 'WDMouse002Sample002',
                                              features = "percent.mt", 
                                              pt.size.factor = 1.8)  + 
                           scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                                              "#2c7fb8","#253494"),
                                                 aesthetics = "fill", 
                                                 limits = c(0.468007,
                                                            9.691789),
                                                 name = "% mit. RNA") +
                           ggtitle('') + theme(legend.position = 'right')
                         ) 

  
  # plotClusters <- SpatialPlot(spatialObjectMergedFiltered,
  #                            features = 'seurat_clusters',
  #                               alpha )  & 
  #   theme(legend.position = "none") & ggtitle('') 
  # 
  # plot
  # 
  # combine HDS on tissue, tissue without spots and box-plot of HDS distribution
  
  layer1plot <- wrap_plots((p1 | p2 | p3 | p4))
  panelplot <- plotJustImages / layer1plot / plotMito 
  panelplot + plot_annotation(
      title = 'HDS Spatial Distribution',
      subtitle = 'Centered score was calculated after merging, filtering and transforming (SCT)\nLiver Cell Atlas spatial transcriptomic data'
    )  
  ggsave('CenteredHDSonTissuePlusTissueSCTFilteredMerged.png',
         path = plotsPath)
  
  # Plot for thesis: 
  
  layer1plot <- wrap_plots(((p1 | p2) / (p3 | p4)))
  
      
  # idea for interpretation, only clusters 0, 3, 4, 7, and 11 are found in 
  # healthy tissue. we can say we associate the other clusters with disease 
  # and look at how the score is distributed in those clusters not present in 
  # healthy tissue 
  
  ## ZOOM on WD samples, by changing range of legend (limits parameter)
  
  min(subset(spatialObjectMergedFiltered, 
             subset = 'orig.ident' == 'WDMouse001Sample001')@meta.data$centeredHDS_SCT)
  
  pp1 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "centeredHDS_SCT",
                           images = 'WDMouse001Sample001',
                           pt.size.factor = 1.8) + 
    scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                       "#2c7fb8","#253494"),
                          aesthetics = "fill", 
                          limits = c(-0.34,0.32),
                          name = "centered HDS") +
    theme(legend.position = 'none'
          # plot.title = element_text(size = 20)
    ) +
    ggtitle("Western Diet\nMouse 1 Sample 1")
  
  pp2 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "centeredHDS_SCT",
                           images = 'WDMouse002Sample001',
                           pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                       "#2c7fb8","#253494"),
                          aesthetics = "fill", 
                          limits = c(-0.34,0.32),
                          name = "centered HDS") +
    theme(legend.position = 'none'
          # plot.title = element_text(size = 20)
    ) +
    ggtitle("Western Diet\nMouse 2 Sample 1")
  
  pp3 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "centeredHDS_SCT",
                           images = 'WDMouse002Sample002',
                           pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  c("#FFFFCC","#a1dab4","#41b6c4",
                                       "#2c7fb8","#253494"),
                          aesthetics = "fill", 
                          limits = c(-0.34,0.32),
                          name = "centered HDS") +
    theme(legend.position = 'right'
          #legend.key.size = unit(1.5, 'cm'),
          #legend.title = element_text(size = 20),
          #plot.title = element_text(size = 20)
    ) +
    ggtitle("Western Diet\nMouse 2 Sample 2")

  
  # IMPORTANT!
  ## pending: use spot annotation info available on liver cell atlas website and
  ## plot on tissue
  
  # beware:not cut  off, ouliers and extreme values included 
  
  plotUmapHDS <- FeaturePlot(spatialObjectMergedFiltered, 
                             'centeredHDS_SCT',
                             pt.size = 0.6, 
                            # min.cutoff = -0.25,
                             #max.cutoff =  0.2, 
                            cols = c("#FFFFCC","#253494")) + 
    theme_dark() 

  plotUmapZones <-DimPlot(spatialObjectMergedFiltered, 
                          reduction = "umap",
                          group.by = 'zonationGroup',
                          cols = c("#a1dab4","#41b6c4",
                                   "#2c7fb8","#253494")
                          )
  
  wrap_plots((umapPlot | umapPlotSplit )/ (plotUmapHDS |plotUmapZones ))
  
  ggsave('UMAPplusCenteredHDSZonationSamples.png',
         path = plotsPath) 
  
  plotUmapMito <- FeaturePlot(spatialObjectMergedFiltered, 
                          'percent.mt_merged',
                          cols = c("#FFFFCC","#253494")
  ) + theme_dark()
  
  plotUmapRbc <- FeaturePlot(spatialObjectMergedFiltered, 
                             'percent.rbc_merged',
                             max.cutoff = 0.1,
                             cols = c("#FFFFCC","#253494") 
  ) + theme_dark()
  
  # cluster 9 has the highest levels of Rbc genes and 
    
    # I noticed something weird. There should not be any values over 1 but there
  # there are values like 6 and so on 
  
  # To Do Friday:
  # plot all of these on tissue
  

  # how can I plot the zones on the tissue ??
  #!!!!!!
  # violin plots of HDS per zone group!! per slice and also grouped 
  
  # First Western Diet Samples All Together
  ZonationHDSViolinPlotWesternDiet <- 
    ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident != "StSpatialSeuratM1S1",], 
           aes(x = zonationGroup, y = centeredHDS)) +
    geom_violin() + xlab("Zones") + ylab("centered HDS") +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12),
          title = element_text(size = 15)
          )  + 
    geom_boxplot(width=0.1, outlier.size = 0) +
    ggtitle('Distribution of Centered HDS per Liver Zone in Western Diet Samples')
  
  
  ggsave(ZonationHDSViolinPlotWesternDiet, 
         path = plotsPath,
         file= 'ZonationHDSViolinPlotWesternDiet.png')
  
  # HDS Violin Plot St Mouse 1 Sample 1
  ZonationHDSViolinPlotM1S1StandardDiet <- 
    ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "StSpatialSeuratM1S1",], 
           aes(x = zonationGroup, y = centeredHDS)) + 
    geom_violin() + xlab("Zones") + ylab("centered HDS") +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12),
          title = element_text(size = 12)
    )  + 
    geom_boxplot(width=0.1, outlier.size = 0) +
    ggtitle('Distribution of Centered HDS per Liver Zone\nin Standard Diet Mouse 1 Sample 1')
  
  # HDS Violin Plot WD Mouse 1 Sample 1
  ZonationHDSViolinPlotM1S1WesternDiet <- 
    ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "WDSpatialSeuratM1S1",], 
           aes(x = zonationGroup, y = centeredHDS)) +
    geom_violin() + xlab("Zones") + ylab("centered HDS") +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12),
          title = element_text(size = 12)
    )  + 
    geom_boxplot(width=0.1, outlier.size = 0) +
    ggtitle('Distribution of Centered HDS per Liver Zone\nin Western Diet Mouse 1 Sample 1')
  
  # HDS Violin Plot WD Mouse 2 Sample 1
  ZonationHDSViolinPlotM2S1WesternDiet <- 
    ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "WDSpatialSeuratM2S1",],
           aes(x = zonationGroup, y = centeredHDS)) +
    geom_violin() + xlab("Zones") + ylab("centered HDS") +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12),
          title = element_text(size = 12)
    )  + 
    geom_boxplot(width=0.1, outlier.size = 0) +
    ggtitle('Distribution of Centered HDS per Liver Zone\nin Western Diet Mouse 2 Sample 1')
  
  # HDS Violin Plot WD Mouse 2 Sample 2
  ZonationHDSViolinPlotM2S2WesternDiet <- 
    ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "WDSpatialSeuratM2S2",], 
           aes(x = zonationGroup, y = centeredHDS)) +
    geom_violin() + xlab("Zones") + ylab("centered HDS") +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12),
          title = element_text(size = 12)
    )  + 
    geom_boxplot(width=0.1, outlier.size = 0) +
    ggtitle('Distribution of Centered HDS per Liver Zone\nin Western Diet Mouse 2 Sample 2')
  
  HDSZonesAllViolins <- (ZonationHDSViolinPlotM1S1StandardDiet|ZonationHDSViolinPlotM1S1WesternDiet)/
    (ZonationHDSViolinPlotM2S1WesternDiet|ZonationHDSViolinPlotM2S2WesternDiet)
  
  ggsave(HDSZonesAllViolins, 
         path = plotsPath,
         file= 'HDSZonesAllViolinPlots.png')
  
  
  
  ########
  # find out what the clusters represent: f.e. which cell types?
  # can I use markers from paper
  #########
  
  ##########
  # Liver Cell Atlas Standard Diet Mice GSE192741
  ##########
  setwd('Desktop/HepatocyteDamageScore/SpatialTranscriptomics/LiverCellAtlasVisium/')
  dataPath <- 'GSE192741_HealthyMiceLivers/'
  plotsPath <- 'GSE192741_HealthyMiceLivers/Plots/'
  
  # GEO data info: 
  
  #All samples are
  # strain: C57Bl/6
  # platform: 10x Visium
  # condition: StSt
  # number of added abs: 0
  # Extracted molecule	total RNA
  
  # Extraction protocol	All Methods listed in Guilliams et al. 
  # Spatial proteogenomics reveals distinct and evolutionarily-conserved 
  # hepatic macrophage niches. Cell. 2022.
  
  # Sample Info
  # 1. Characteristics	shortfilename: 
  # JBO001
  # number of spots: 1646
  # 2. Characteristics	shortfilename: JBO002
  # number of spots: 1651
  # 3. Characteristics	shortfilename: 
  # number of spots: 
  # 4. Characteristics	shortfilename: 
  # number of spots: 
  
  
  
  annotationNAFLD <- read.csv(file = paste0(dataPath,"annot_mouseNafldVisium.csv"))
  
  
  spatialSamples <- list(
    Load10X_Spatial(data.dir = paste0(dataPath,'StSt/Mouse001Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5",
                    assay = "LiverCellAtlasSpatial",
                    slice = "StMouse001Sample001",
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "StSt/Mouse001Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse001Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5", 
                    assay = "LiverCellAtlasSpatial", 
                    slice = "WDMouse001Sample001", 
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse001Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse002Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5",
                    assay = "LiverCellAtlasSpatial",
                    slice = "WDMouse002Sample001",
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse002Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse002Sample002/'),
                    filename = "filtered_feature_bc_matrix.h5", 
                    assay = "LiverCellAtlasSpatial", 
                    slice = "WDMouse002Sample002", 
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse002Sample002/")))
  )
  
  sampleNames <- c("StSpatialSeuratM1S1",
                   "WDSpatialSeuratM1S1",
                   "WDSpatialSeuratM2S1", 
                   "WDSpatialSeuratM2S2")
  
  names(spatialSamples) <- sampleNames
  
  sampleCoordinates <- list(
    spatialSamples[[1]]@images$StMouse001Sample001@coordinates,
    spatialSamples[[2]]@images$WDMouse001Sample001@coordinates,
    spatialSamples[[3]]@images$WDMouse002Sample001@coordinates,
    spatialSamples[[4]]@images$WDMouse002Sample002@coordinates
  )
  
  
}
