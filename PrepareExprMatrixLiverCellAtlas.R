# Expression Matrices to Apply HDS # 
# # run on cluster # #

# Liver Cell Atlas:
# keeping only hepatocytes and single nuclei RNAseq


library(Seurat)
library(sctransform)
library(ggplot2)
library(tidyverse)

# First Standard Diet Counts

pathData = 
  '~/.../rawData_mouseStSt/countTable_mouseStSt/'
pathAnnotation = 
  '~/.../annot_mouseStStAll.csv'

outputPath = '~/...'

annotation.mouseStSt <- read.csv(
  file = pathAnnotation)

# read raw counts 
counts <- ReadMtx(mtx =  paste0(pathData,'matrix.mtx.gz'), 
                  cells = paste0(pathData,'barcodes.tsv.gz'),
                  features = paste0(pathData,'features.tsv.gz'),
                  feature.column = 1) 

# create seurat object 
lca.mouse.hep.nucSEq.stst <- CreateSeuratObject(counts = counts, 
                                                project = "LCAStdDiet",
                                                min.cells = 3)

remove(counts)

# save indices to keep only single nuclei hepatocyte reads
cellsToKeep <- annotation.mouseStSt[annotation.mouseStSt$annot == 'Hepatocytes' & 
                                      annotation.mouseStSt$typeSample == 'nucSeq', 
                                    'cell'] 
# subset seurat object
lca.mouse.hep.nucSEq.stst <- subset(lca.mouse.hep.nucSEq.stst, 
                                    cells = cellsToKeep)

remove(cellsToKeep)

# add selected metadata to seurat object
annotation.mouseStSt <- annotation.mouseStSt[annotation.mouseStSt$annot == 'Hepatocytes' & 
                                               annotation.mouseStSt$typeSample == 'nucSeq', ] 
annotation.mouseStSt <- annotation.mouseStSt[is.element(annotation.mouseStSt$cell, 
                                                        names(lca.mouse.hep.nucSEq.stst$orig.ident)),] 
lca.mouse.hep.nucSEq.stst@meta.data <- annotation.mouseStSt


# repeat steps for NAFLD counts
pathData = 
  '~/.../LiverCellAtlas/rawData_mouseNafld/countTable_mouseNafld/'
pathAnnotation = 
  '~/.../annot_mouseNafldAll.csv'

annotation.mouseNAFLD <- read.csv(
  file = pathAnnotation)

# read raw counts 
counts <- ReadMtx(mtx =  paste0(pathData,'matrix.mtx.gz'), 
                  cells = paste0(pathData,'barcodes.tsv.gz'),
                  features = paste0(pathData,'features.tsv.gz'),
                  feature.column = 1) 

# create seurat object 
lca.mouse.hep.nucSEq.nafld <- CreateSeuratObject(counts = counts, 
                                                 project = "LCANafld",
                                                 min.cells = 3)

remove(counts)

# save indices to keep only single nuclei hepatocyte reads
cellsToKeep <- annotation.mouseNAFLD[annotation.mouseNAFLD$annot == 'Hepatocytes' &
                                       annotation.mouseNAFLD$typeSample == 'nucSeq', 
                                     'cell'] 
# subset seurat object
lca.mouse.hep.nucSEq.nafld <- subset(lca.mouse.hep.nucSEq.nafld, 
                                     cells = cellsToKeep)
remove(cellsToKeep)

# add selected meta data to seurat object 

annotation.mouseNAFLD <- annotation.mouseNAFLD[annotation.mouseNAFLD$annot == 'Hepatocytes' & 
                                                 annotation.mouseNAFLD$typeSample == 'nucSeq', ] 
annotation.mouseNAFLD <- annotation.mouseNAFLD[is.element(annotation.mouseNAFLD$cell, 
                                                          names(lca.mouse.hep.nucSEq.nafld$orig.ident)),] 
lca.mouse.hep.nucSEq.nafld@meta.data <- annotation.mouseNAFLD

# merge Sdt. Diet and NAFLD counts 

merged.mouse.lca.hep.nucSeq <-
  merge(x = temp_stst,
        y = temp_nafld, 
        add.cell.ids = c("StSt", "Nafld"), 
        project = "LiverCellAtlas")

# short QC and filtering


merged.mouse.lca.hep.nucSeq[["percent.mt"]] <- PercentageFeatureSet(merged.mouse.lca.hep.nucSeq, 
                                                                    pattern = "^mt-")


# Filtering
merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, 
                                      subset = nFeature_RNA > 500 &
                                        nFeature_RNA < 6000 &
                                        percent.mt <= 5 & 
                                        nCount_RNA < 30000 )

# Save expression matrix to apply HDS on
exprMatrix <- GetAssayData(object = merged.mouse.lca.hep.nucSeq , 
                             slot = 'counts')

temp <- merged.mouse.lca.hep.nucSeq@meta.data
annotations <- data.frame('Sample' = temp$sample , 
                          'cell_id' = rownames(temp),
                          'condition' = gsub(pattern = "_.*", 
                                             replacement = "", 
                                             rownames(temp)) )

exprMatrix <- as.matrix(exprMatrix)
saveRDS(exprMatrix, 
        'ExprMatrixMouseMergedHepatocytesnucSeqQCfiltered.rds')
write.csv(annotations, 
          'ExprMatrixMouseMergedHepatocytesnucSeqQCfilteredAnnotations.csv')

remove(exprMatrix, merged.mouse.lca.hep.nucSeq, annotations, temp,
       annotation.mouseNAFLD, annotation.mouseStSt)

#### dont need anything above from this comment really

# Save expression matrix to apply HDS on
merged.mouse.lca.hep.nucSeq <- readRDS('scRNAseq/Liver Cell Atlas Data/mouseMergedHepatocytesnucSeqQCfitleredSCT.rds')

exprMatrix <- GetAssayData(object = merged.mouse.lca.hep.nucSeq , 
                           slot = 'counts')

temp <- merged.mouse.lca.hep.nucSeq@meta.data
annotations <- data.frame('Sample' = temp$sample , 
                          'cell_id' = rownames(temp),
                          'condition' = gsub(pattern = "_.*", 
                                             replacement = "", 
                                             rownames(temp)) )

exprMatrix <- as.matrix(exprMatrix)
saveRDS(exprMatrix, 
        'scRNAseq/Liver Cell Atlas Data/ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv')
write.csv(annotations, 
          'scRNAseq/Liver Cell Atlas Data/ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv')

remove(exprMatrix, merged.mouse.lca.hep.nucSeq, annotations, temp,
       annotation.mouseNAFLD, annotation.mouseStSt)


