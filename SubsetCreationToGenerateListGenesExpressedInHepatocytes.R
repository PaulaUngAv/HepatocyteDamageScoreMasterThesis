library(Seurat)

## merge

temp_stst <- readRDS(file = 'countTable_mouseStSt_Hepatocytes_nucSeq.rds')
annotation_temp <- read.csv(file ='annot_mouseStStAll.csv')
annotation_temp <- annotation_temp[annotation_temp$annot == 'Hepatocytes', ]
annotation_temp<- annotation_temp[is.element(annotation_temp$cell, names(temp_stst$orig.ident)),] 
temp_stst@meta.data$sample <- annotation_temp$sample
temp_stst@meta.data$digest <- annotation_temp$digest
temp_stst@meta.data$typeSample <- annotation_temp$typeSample
temp_stst@meta.data$celltype<- annotation_temp$annot

temp_nafld <- readRDS(file = 'countTable_mouseNafld_Hepatocytes_nucSeq.rds')
annotation_temp <- read.csv(file ='annot_mouseNafldAll.csv')
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