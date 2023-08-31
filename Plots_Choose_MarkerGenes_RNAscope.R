
## Load Data 
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)

{
  path <- '~/Desktop/Project Module Beyer /'
  studies <- c('GSE99010/GSE99010_CCl4_deg.csv',
             'GSE99010/GSE99010_WesternDiet_deg.csv',
             'GSE97234/GSE97234_AH_deg.csv',
             'GSE97234/GSE97234_ASH_deg.csv',
             'GSE97024/GSE97024_Oxr_BrainKO_deg.csv',
             'GSE83240/GSE83240_DEN_deg.csv',
             'GSE61117/GSE61117_Trim24_KO_deg.csv',
             'GSE153580/GSE153580_STZ_deg.csv',
             'GSE148849/GSE148849_fastfooddiet_deg.csv',
             'GSE138419/GSE138419_AMLN_deg.csv',
             'GSE137449/GSE137449_CDAHFD_diet_deg.csv',
             'GSE135050/GSE135050_hfcfdiet_deg.csv',
             'GSE132040/GSE132040_young(6_9_12)vs_old(18_21_24)_deg.csv',
             'GSE124694/GSE124694_Arid1a_LiverKO_deg.csv',
             'GSE123894/GSE123894_Fructose_only_DBA2Jstrain_deg.csv',
             'GSE119953/GSE119953_HFD_deg.csv',
             'GSE119953/GSE119953_DDC_deg.csv',
             'GSE119441/GSE119441_HFD_deg.csv',
             'GSE119441/GSE119441_PFOA_deg.csv',
             'GSE114261/GSE114261_STAM_deg.csv',
             'GSE111828/GSE111828_Acetaminophen_deg.csv',
             'GSE109431/GSE109431_LHR1_KO_deg.csv',
             'GSE106369/GSE106369_MOFKO_deg.csv'
  )
  
  DElist <- lapply( seq( studies ), function(ii){
    read.csv(paste0(path,studies[[ii]]))
  })

  names(DElist) <- sub( "_deg.csv","",basename(studies))
  UniHepDamageSignature <- read.csv(file = '~/Desktop/Project Module Beyer /Universal_Hepatocyte_Damage_Signature_9June22.csv')
  
  LFC_perModel_base <- data.frame('Gene_Symbol' = head(UniHepDamageSignature$gene_symbol, 42),
                             'Average_LFC' = head(UniHepDamageSignature$direction_foldchange, 42))
  
  LFCs <- lapply(seq(names(DElist)), function(ii){
    
    temp <- DElist[[ii]]
    LFC <- temp[temp$Gene.Symbol %in% LFC_perModel_base$Gene_Symbol, c('Gene.Symbol','log2FoldChange')]
    LFC <- as.vector(LFC[order(LFC$Gene.Symbol),'log2FoldChange'])
    ## BUG HERE !! I think 
    return(LFC)
  })
  
  names(LFCs) <- names(DElist)
  LFC_perModel <- cbind(LFC_perModel_base, do.call(cbind, LFCs))
  LFC_perModel_matrix <- LFC_perModel 
  rownames(LFC_perModel_matrix) <- LFC_perModel$Gene_Symbol
  LFC_perModel_matrix <- as.matrix(LFC_perModel_matrix[,3:25])
  
  # heatmap base R
  heatmap(LFC_perModel_matrix,scale = 'none')
  library(circlize)
  col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
  col_fun(seq(-1, 1))
  ordered_cols <- c("GSE97234_AH","GSE97234_ASH","GSE99010_CCl4","GSE119953_DDC","GSE153580_STZ",
                    "GSE119441_PFOA","GSE111828_Acetaminophen","GSE119441_HFD",
                    "GSE119953_HFD", "GSE148849_fastfooddiet",
                    "GSE138419_AMLN",
                    "GSE99010_WesternDiet","GSE135050_hfcfdiet","GSE137449_CDAHFD_diet",
                    "GSE83240_DEN","GSE114261_STAM",
                    "GSE61117_Trim24_KO","GSE97024_Oxr_BrainKO" ,"GSE106369_MOFKO",
                    "GSE124694_Arid1a_LiverKO","GSE109431_LHR1_KO", "GSE132040_young(6_9_12)vs_old(18_21_24)",
                    "GSE123894_Fructose_only_DBA2Jstrain")
  col1 <- HeatmapAnnotation(model = c('CCl4','HFCF', 'EtOH','EtOH', 'KO' ,'DEN',
                                      'KO', 'STZ','HFCF', 'HFCF', 'CDAHFD','HFCF',
                                      'Aging','KO', 'Fructose','HFD', 'DDC', 'HFD',
                                      'PFOA', 'STAM', 'Acetaminophen','KO','KO'),
                            which = 'col'
                            )
  Heatmap(LFC_perModel_matrix, name = 'LFC' , col = col_fun,
          column_title = "Disease model and GEO identifier",
          column_title_side = "top",
          row_title = "Genes",
          column_order = ordered_cols,
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          top_annotation = col1,
          show_column_names = T
          )
  
  Heatmap(LFC_perModel_matrix, name = 'LFC' , col = col_fun,
          column_title = "Disease model and GEO identifier",
          column_title_side = "top",
          row_title = "Genes",
          cluster_columns = T,
          cluster_rows = T,
          show_column_names = T

          
  )
  

##### will not use this 
  ## log fold change table 
  LFCs_binary <- lapply(seq(LFCs), function(ii){
    return(sign(LFCs[[ii]]))
  })
  names(LFCs_binary) <- names(DElist)
  LFC_perModel_binary <- cbind(LFC_perModel_base, do.call(cbind, LFCs_binary))
  LFC_perModel_binary_forplot <- pivot_longer(
    data = LFC_perModel_binary,
    cols = !c(Gene_Symbol,Average_LFC),
    names_to = "Dataset",
    values_to = "LFC"
    
  )
  
  LFC_perModel_binary_matrix <- LFC_perModel_binary
  rownames(LFC_perModel_binary_matrix) <- LFC_perModel_binary$Gene_Symbol
  LFC_perModel_binary_matrix <- as.matrix(LFC_perModel_binary_matrix[,-c(1,2)])
  
  #Binary Heatmap
  Heatmap(LFC_perModel_binary_matrix, name = 'LFC' , col = col_fun,
          column_title = "Disease model and GEO identifier",
          column_title_side = "top",
          row_title = "Genes",
          cluster_columns = T,
          cluster_rows = T,
          show_column_names = T
          
          
  )

  
  ## KEGG 
  

  library(gage)
  library(pathview)
  library(gageData)
  
  data(kegg.sets.mm)
  pathview()
  
  ## gene expression 
  library(Seurat)
  merged_liver_cell_atlas_hepnuc <- readRDS('../scRNAseq/Liver Cell Atlas Data/countTable_mouseMerged_Hepatocytes_nucSeq.rds')
  DotPlot(merged_liver_cell_atlas_hepnuc, features = LFC_perModel$Gene_Symbol) + RotatedAxis()
  remove(merged_liver_cell_atlas_hepnuc)
  
  }
