merged.mouse.lca.hep.nucSeq <- readRDS('scRNAseq/Liver Cell Atlas Data/countTable_mouseMerged_Hepatocytes_nucSeq.rds')
VlnPlot(merged.mouse.lca.hep.nucSeq, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt" ), 
        ncol = 3,pt.size = F) 


VlnPlot(merged.mouse.lca.hep.nucSeq,
        features = c("percent.mt" ), pt.size = F, group.by = 'typeSample' ) 

plot1 <- FeatureScatter(merged.mouse.lca.hep.nucSeq, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot1
merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, 
                                      subset = nFeature_RNA > 200 &
                                        nFeature_RNA < 6000 &
                                        percent.mt <= 5 & 
                                        nCount_RNA < 30000 )
exprMatrices <- GetAssayData(object = merged.mouse.lca.hep.nucSeq , slot = 'counts')
temp <- merged.mouse.lca.hep.nucSeq@meta.data
annotations <- data.frame('Sample' = temp$sample , 
                          'cell_id' = rownames(temp),
                          'condition' = gsub(pattern = "_.*", replacement = "", rownames(temp)) )

exprMatrices <- as.matrix(exprMatrices)
saveRDS(exprMatrices, 'scRNAseq/Liver Cell Atlas Data/ExprMatrix_QCfiltered_mouseMerged_Hepatocytes_nucSeq.rds')
write.csv(annotations, 'scRNAseq/Liver Cell Atlas Data/ExprMatrix_QCfiltered_mouseMerged_Hepatocytes_nucSeq_annotations.csv')
