## Gene Expressed in Hepatocytes ##
# final version 

# 1, Subset Liver Cell Atlas Counts to only Hepatocyte Single Nuclei Counts
library(Seurat)
library(sctransform)
library(ggplot2)
library(tidyverse)

# First Standard Diet Counts

pathData = 
   '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/'
pathAnnotation = 
  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseStStAll.csv'

outputPath = '/data/public/pungerav/liver_disease_score/data_disease_score/output/GenesExpressedInHepatocytes/'

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
print(cellsToKeep)

remove(cellsToKeep)

# add selected metadata to seurat object
annotation.mouseStSt <- annotation.mouseStSt[annotation.mouseStSt$annot == 'Hepatocytes' &
                                               annotation.mouseStSt$typeSample == 'nucSeq', ]
annotation.mouseStSt <- annotation.mouseStSt[is.element(annotation.mouseStSt$cell,
                                                        names(lca.mouse.hep.nucSEq.stst$orig.ident)),]
lca.mouse.hep.nucSEq.stst@meta.data <- lca.mouse.hep.nucSEq.stst@meta.data[annotation.mouseStSt$cell,]

if (identical(rownames(lca.mouse.hep.nucSEq.stst@meta.data),
              annotation.mouseStSt$cell) == TRUE){
  lca.mouse.hep.nucSEq.stst@meta.data <- cbind(lca.mouse.hep.nucSEq.stst@meta.data,
                                               annotation.mouseStSt)
}

# repeat steps for NAFLD counts
pathData = 
  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/'
pathAnnotation = 
  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseNafldAll.csv'

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
print(cellsToKeep)
# subset seurat object
lca.mouse.hep.nucSEq.nafld <- subset(lca.mouse.hep.nucSEq.nafld, 
                                    cells = cellsToKeep)
remove(cellsToKeep)

# add selected metadata to seurat object
annotation.mouseNAFLD <- annotation.mouseNAFLD[annotation.mouseNAFLD$annot == 'Hepatocytes' &
                                               annotation.mouseNAFLD$typeSample == 'nucSeq', ]
annotation.mouseNAFLD <- annotation.mouseNAFLD[is.element(annotation.mouseNAFLD$cell,
                                                        names(lca.mouse.hep.nucSEq.nafld$orig.ident)),]
lca.mouse.hep.nucSEq.nafld@meta.data <- lca.mouse.hep.nucSEq.nafld@meta.data[annotation.mouseNAFLD$cell,]

if (identical(rownames(lca.mouse.hep.nucSEq.nafld@meta.data),
              annotation.mouseNAFLD$cell) == TRUE){
  lca.mouse.hep.nucSEq.nafld@meta.data <- cbind(lca.mouse.hep.nucSEq.nafld@meta.data,
                                               annotation.mouseNAFLD)
}

# merge Sdt. Diet and NAFLD counts 

merged.mouse.lca.hep.nucSeq <-
  merge(x = lca.mouse.hep.nucSEq.stst,
        y = lca.mouse.hep.nucSEq.nafld, 
        add.cell.ids = c("StSt", "Nafld"), 
        project = "LiverCellAtlas")

# short QC and filtering


merged.mouse.lca.hep.nucSeq[["percent.mt"]] <- PercentageFeatureSet(merged.mouse.lca.hep.nucSeq, 
                                                                    pattern = "^mt-")

DimsBeforeFiltering <- dim(merged.mouse.lca.hep.nucSeq)
print(DimsBeforeFiltering)

# Plots pre Filtering
prefilterFeauturePlot <- VlnPlot(merged.mouse.lca.hep.nucSeq, 
                                 features = c("nFeature_RNA", 
                                              "nCount_RNA",
                                              "percent.mt" ), 
                                 ncol = 3,pt.size = FALSE) 

# Filtering
merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, 
                                      subset = nFeature_RNA > 500 &
                                        nFeature_RNA < 6000 &
                                        percent.mt <= 2.5 & 
                                        nCount_RNA < 30000 )

# post-filtering plots
prefilterFeauturePlot <- VlnPlot(merged.mouse.lca.hep.nucSeq, 
                                 features = c("nFeature_RNA", 
                                              "nCount_RNA",
                                              "percent.mt" ), 
                                 ncol = 3,pt.size = FALSE) 

# SCT transform counts

merged.mouse.lca.hep.nucSeq <- SCTransform(merged.mouse.lca.hep.nucSeq,
                                           verbose = FALSE)



## Select Genes


# Detection Rate in hepatocytes (LFC between single nucleus and single cell)
foldchange <- FoldChange(merged.mouse.lca.hep.nucSeq,
                         ident.1 = 'nucSeq', 
                         group.by = 'typeSample',
                         assay = "SCT")

detection_rate_df <- data.frame('Genes' = row.names(foldchange), 
                                'detection_rate' = foldchange[,2],
                                'log2_dec_rate' = log2(foldchange[,2]+ 0.0001))

ggplot(data= detection_rate_df, aes(x=detection_rate)) + geom_histogram(binwidth = 0.005)
ggplot(data = detection_rate_df, aes(x = log2_dec_rate)) + geom_density()

summaryStats <- summary(detection_rate_df$log2_dec_rate, na.rm = TRUE)
quantilesLog2DetRate <- quantile(detection_rate_df$log2_dec_rate)

# Plotting theshold 
plotDetectionRate <- ggplot(data= detection_rate_df, aes(x=log2_dec_rate)) +
  geom_histogram(binwidth = 0.05) + 
  geom_vline(xintercept = c(quantilesLog2DetRate[2], 
                            quantilesLog2DetRate[3],
                            quantilesLog2DetRate[4]),
             colour = 'red') + 
  xlab('log2(detection rate of gene + 0.0001)') + 
  ylab('number of genes') + 
  geom_text(aes(x=quantilesLog2DetRate[2], 
                label="1. Quartile", 
                y=1100), 
            colour="red", 
            angle=90, 
            vjust = -1) +
  geom_text(aes(x=quantilesLog2DetRate[3], 
                label="2. Quartile", 
                y=1100), 
            colour="red", 
            angle=90, 
            vjust = -1) +
  geom_text(aes(x=quantilesLog2DetRate[4], 
                label="3. Quartile", 
                y=1100), 
            colour="red", 
            angle=90,
            vjust = -1)

# Save Output

# save gene list 

detection_rate_df <- 
  detection_rate_df[detection_rate_df$log2_dec_rate >  -5, ]

detection_rate_df <- 
  detection_rate_df[order(detection_rate_df$log2_dec_rate, decreasing = TRUE),]

rownames(detection_rate_df) <- 1:length(detection_rate_df$Genes)

# Save filtered filtered Seurat object
# Seurat object
saveRDS(merged.mouse.lca.hep.nucSeq, 
        paste0(outputPath,'mouseMergedHepatocytesnucSeqQCfitleredSCT.rds'))



# Save gene list

write.csv(detection_rate_df, 
          file = paste0(outputPath, 'GenesExpressedInHepatocytes.csv'))

# Save plots

ggsave(
  'prefilterFeauturePlot.png',
  plot = prefilterFeauturePlot,
  path = outputPath
)

ggsave(
  'prefilterFeauturePlot.png',
  plot = prefilterFeauturePlot,
  path = outputPath
)

ggsave(
  'plotDetectionRate.png',
  plot = plotDetectionRate,
  path = outputPath
)

# Save session info
write.table(toLatex(sessionInfo()), file = 'sessionInfo.txt', sep = '/')






