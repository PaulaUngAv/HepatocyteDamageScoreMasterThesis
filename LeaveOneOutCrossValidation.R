#### ==================#####
####  Cross-Validation  #### 
#### ==================#####

# Load packages and data
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
  library(patchwork)
  
  #### load functions
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

  # loads list of all DEG tables
  DElist <- lapply( seq(studies), function(ii){
    read.csv(paste0(path,"AllBulkStudiesDEGAnalysis/",studies[[ii]]))}
  )
  names(DElist) <- basename(studies)
  
  HUDS <- 
    read.csv(file = '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')
  
  geneFilter <- 
    read.csv(file ='~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')
  geneHepa <- geneFilter$Genes
  
  ### load bulk expression matrices
  ll <- list.files( full.names = T, path= 
                      "~/Desktop/HepatocyteDamageScore/All_Count_Tables")
  ### load annotations of bulk expression matrices
  aa <- list.files(full.names = T, path = 
                     "~/Desktop/HepatocyteDamageScore/All_Annotations")
  bulk_explist <- lapply( ll, FUN = read.delim, sep = ' ')
  bulk_explist <- lapply(seq(bulk_explist), function(ii){
    if (length(unique(bulk_explist[[ii]][,1])) != nrow(bulk_explist[[ii]])){
      print(ll[ii])
      print(ii)
    }else{
      rownames(bulk_explist[[ii]]) <- bulk_explist[[ii]][,1]
      bulk_explist[[ii]] <- data.matrix(bulk_explist[[ii]][,-1])
    }
  })
  bulk_annlist <- lapply( aa, FUN = read.csv, header = T) 
  
  ### reformat names of data frames in lists to only contain GEO IDs
  names(bulk_explist) <-sub( "_.*","",basename(ll))
  names(bulk_annlist) <- sub("SraRunTable_","",sub( ".txt","",basename(aa)))
  
  ### Creates (temporal) list of vectors. 
  # Vector contains sample IDs corresponding to study GEO ID
  annotation_sample_ids <- lapply(seq(bulk_explist), function(ii){
    colnames(bulk_explist[[ii]])
  })
  names(annotation_sample_ids) <- names(bulk_explist)
 
   ### Creates list of data frames that include the condition of each samples 
  ## per study (GEO identifier)
  condition_sample_list <- list(data.frame('sample_id' = annotation_sample_ids[[1]],
                                           'condition' =  bulk_annlist[[1]]$Genotype[bulk_annlist[[1]]$Run ==  annotation_sample_ids[[1]]]),
                                data.frame('sample_id' =annotation_sample_ids[[2]],
                                           'condition' = bulk_annlist[[2]]$genotype.variation[bulk_annlist[[2]]$Run ==  annotation_sample_ids[[2]]]),
                                data.frame('sample_id' = annotation_sample_ids[[3]],
                                           'condition' = bulk_annlist[[3]]$time[bulk_annlist[[3]]$Run ==  annotation_sample_ids[[3]]]),
                                data.frame('sample_id' = annotation_sample_ids[[4]],
                                           'condition' = bulk_annlist[[4]]$mouse_model[bulk_annlist[[4]]$Run  %in% annotation_sample_ids[[4]]]),
                                data.frame('sample_id' = annotation_sample_ids[[5]],
                                           'condition' = bulk_annlist[[5]]$treatment[bulk_annlist[[5]]$Run ==  annotation_sample_ids[[5]]]),
                                data.frame('sample_id' = annotation_sample_ids[[6]],
                                           'condition' = paste(bulk_annlist[[6]]$treatment_1[bulk_annlist[[6]]$Sample.Name ==  annotation_sample_ids[[6]]],
                                                               bulk_annlist[[6]]$treatment_2[bulk_annlist[[6]]$Sample.Name ==  annotation_sample_ids[[6]]])),
                                data.frame('sample_id' = annotation_sample_ids[[7]],
                                           'condition' = bulk_annlist[[7]]$treatment[bulk_annlist[[7]]$Run ==  annotation_sample_ids[[7]]]),
                                data.frame('sample_id' = annotation_sample_ids[[8]],
                                           'condition' = bulk_annlist[[8]]$genotype.variation[bulk_annlist[[8]]$Run ==  annotation_sample_ids[[8]]]),
                                data.frame('sample_id' = annotation_sample_ids[[9]],
                                           'condition' = bulk_annlist[[9]]$characteristics..age[bulk_annlist[[9]]$Sample.name %in% annotation_sample_ids[[9]]] ),
                                data.frame('sample_id' = annotation_sample_ids[[10]],
                                           'condition' = paste(bulk_annlist[[10]]$diet[bulk_annlist[[10]]$Run ==  annotation_sample_ids[[10]]],
                                                               bulk_annlist[[10]]$months[bulk_annlist[[10]]$Run ==  annotation_sample_ids[[10]]],
                                                               "month", sep = '_')),
                                data.frame('sample_id' = annotation_sample_ids[[11]],
                                           'condition' = paste(bulk_annlist[[11]]$diet, bulk_annlist[[11]]$sex)),
                                
                                data.frame('sample_id' = annotation_sample_ids[[12]],
                                           'condition' = paste(bulk_annlist[[12]]$Genotype[bulk_annlist[[12]]$Run ==  annotation_sample_ids[[12]]],
                                                               bulk_annlist[[12]]$diet[bulk_annlist[[12]]$Run ==  annotation_sample_ids[[12]]], sep = " ")),
                                data.frame('sample_id' = annotation_sample_ids[[13]],
                                           'condition' = paste(bulk_annlist[[13]][1:18,]$diet,bulk_annlist[[13]][1:18,]$treatment )),
                                data.frame('sample_id' = annotation_sample_ids[[14]],
                                           'condition' = c(paste(bulk_annlist[[14]]$disease_state[1:9],
                                                                 bulk_annlist[[14]]$time[1:9],sep = ' '), 
                                                           paste(rep('STZ', 9), bulk_annlist[[14]]$time[10:18]))),
                                data.frame('sample_id' = annotation_sample_ids[[15]],
                                           'condition' = bulk_annlist[[15]]$Genotype),
                                data.frame('sample_id' = annotation_sample_ids[[16]],
                                           'condition' = bulk_annlist[[16]]$source_name[bulk_annlist[[16]]$Run %in% annotation_sample_ids[[16]]]),
                                data.frame('sample_id' = annotation_sample_ids[[17]],
                                           'condition' = bulk_annlist[[17]]$Genotype),
                                data.frame('sample_id' = annotation_sample_ids[[18]],
                                           'condition' = bulk_annlist[[18]]$disease_state),
                                data.frame('sample_id' = annotation_sample_ids[[19]],
                                           'condition' = paste(bulk_annlist[[19]]$ccl4_treatment,
                                                               bulk_annlist[[19]]$diet,bulk_annlist[[19]]$diet_duration,
                                                               bulk_annlist[[19]]$tissue))
  )
  
  names(condition_sample_list) <- names(bulk_annlist)

  
  
  DEGbyModel <- list("EtOH" = c('GSE97234_AH_deg.csv',
                                'GSE97234_ASH_deg.csv'),
                     "CCl4" = 'GSE99010_CCl4_deg.csv',
                     "DDC" =  'GSE119953_DDC_deg.csv',
                     "STZ" =  'GSE153580_STZ_deg.csv',
                     "PFOA" = 'GSE119441_PFOA_deg.csv',
                     "DEN" = 'GSE83240_DEN_deg.csv',
                     "Acetaminophen" = 'GSE111828_Acetaminophen_deg.csv',
                     "High Fat Diet" = c( 'GSE119953_HFD_deg.csv',
                                          'GSE119441_HFD_deg.csv'),
                     "Fructose" = 'GSE123894_Fructose_only_DBA2Jstrain_deg.csv',
                     "HFCF" = c('GSE135050_hfcfdiet_deg.csv',
                                'GSE99010_WesternDiet_deg.csv',
                                'GSE148849_fastfooddiet_deg.csv',
                                'GSE138419_AMLN_deg.csv'),
                     "CDAD-HFD" = 'GSE137449_CDAHFD_diet_deg.csv',
                     "STAM" = 'GSE114261_STAM_deg.csv',
                     "Aging" = 'GSE132040_young(6_9_12)vs_old(18_21_24)_deg.csv'
  )
  
  GEOidentifierByModel <- list("EtOH" = 'GSE97234',
                               "CCl4" = 'GSE99010',
                               "DDC" =  'GSE119953',
                               "STZ" =  'GSE153580',
                               "PFOA" = 'GSE119441',
                               "DEN" = 'GSE83240',
                               "Acetaminophen" = 'GSE111828',
                               "High Fat Diet" = c('GSE119953',
                                                   'GSE119441'),
                               "Fructose" = 'GSE123894',
                               "HFCF" = c('GSE135050',
                                          'GSE99010',
                                          'GSE148849',
                                          'GSE138419'),
                               "CDAD-HFD" = 'GSE137449',
                               "STAM" = 'GSE114261',
                               "Aging" = 'GSE132040'
  )  
  
}

#### cross-validate leave one out 
{
  
  # leave leave one model out
  HDS_cv.lomo <- lapply(names(DEGbyModel), function(model){
    print(model)
    print(!names(DElist)%in%DEGbyModel[[model]])
    DSig_cv.onemodel <- damage_signature.func( DElist[ !names(DElist)%in%DEGbyModel[[model]]], 
                                               geneFilter = geneHepa)
    HDS_cv.withoutmodel <- lapply( GEOidentifierByModel[[model]], function(identifier){
      print(identifier)
      # skip names if expression and DE names don't match
      if( identifier %in% names(bulk_explist) ) {
        print(identifier %in% names(bulk_explist))
        DS_calc.func( exprMatrices = bulk_explist[[identifier]],
                      DSignature= DSig_cv.onemodel, 
                      ntop=42 )
      }
    })
    
  })
}


 ## !!  calculate HDS for all of them without taking any study out to compare

# remove empty elements
HDS_cv.lomo <- HDS_cv.lomo[!sapply(HDS_cv.lomo, is.null)]
names( HDS_cv.lomo ) <- names(GEOidentifierByModel)

# preparing results for ggplot

HDS_cv.lomo.ListDataFrames <- lapply( seq(HDS_cv.lomo) , function(ii){
  
  if(length(HDS_cv.lomo[[ii]]) == 1) {
    model <- names(HDS_cv.lomo)[ii]
    temp_id <- GEOidentifierByModel[[model]]
    print(model)
    print(temp_id)
    datt <- cbind.data.frame( score = unname(HDS_cv.lomo[[ii]][[1]]) ,
                              model = rep(model, length(HDS_cv.lomo[[ii]][[1]])) ,
                              condition = as.factor( condition_sample_list[[temp_id]]$condition[ match( names(HDS_cv.lomo[[ii]][[1]]),
                                                                                                        condition_sample_list[[temp_id]]$sample_id)]
                              ))
  } else {
    model <- names(HDS_cv.lomo)[ii]
    temp_id <- GEOidentifierByModel[[model]]
    print(model)
    print(temp_id)
    list_datt <- lapply(seq(GEOidentifierByModel[[model]]), function(jj){
      datt <- cbind.data.frame( score = unname(HDS_cv.lomo[[ii]][[jj]]) ,
                                model = rep(model, length(HDS_cv.lomo[[ii]][[jj]])) ,
                                condition = as.factor( condition_sample_list[[temp_id[jj]]]$condition[ match( names(HDS_cv.lomo[[ii]][[jj]]),
                                                                                                              condition_sample_list[[temp_id[jj]]]$sample_id)]
                                ))
    })
    datt <- rbind.fill(list_datt)
    
  }
  #datt$condition <- relevel( datt$condition , ref = "StSt")
  
  return(datt)
}
)

# first 6 models of GeoIndetifier and model 9 (STAM)
FirstGroupHDSListDataFrames <- lapply(c(seq(1,6),12), function(ii){
  temp <- HDS_cv.lomo.ListDataFrames[[ii]]
  temp$GEOid <- rep(GEOidentifierByModel[[ii]][1],
                    length(HDS_cv.lomo.ListDataFrames[[ii]]$score))
  return(temp)
  
})

# Western Diet, High-Fat, High-Fructose, CDAD-HFD
SecondGroupHDS <- lapply(seq(8,11), function(ii){
  temp <- HDS_cv.lomo.ListDataFrames[[ii]]
  if(length(GEOidentifierByModel[[ii]]) > 1 ){
    # this will be changed then manually
    temp$GEOid <- rep(toString(GEOidentifierByModel[[ii]]), 
                      length(HDS_cv.lomo.ListDataFrames[[ii]]$score))
  }else{
    temp$GEOid <- 
      rep(GEOidentifierByModel[[ii]][1], 
          length(HDS_cv.lomo.ListDataFrames[[ii]]$score))
    }
  
  return(temp)
  
}
)

outputPath <- '~/Desktop/HepatocyteDamageScore/Results/CrossValidation/'
write.csv(Reduce(FirstGroupHDSListDataFrames,f = rbind), 
          file = paste0(outputPath, 'LeaveOneOutCrossValidationFirstGroup.csv'))
write.csv(Reduce(SecondGroupHDS,f = rbind), 
          file = paste0(outputPath, 'LeaveOneOutCrossValidationSecondGroup.csv'))

## edited all prior tables on excel to manually label different conditions to 
## groups of controls and treatment for plotting 
FirstGroupHDSListDataFrames <- 
  read.csv(file = paste0(outputPath,
                         'LeaveOneOutCrossValidationFirstGroupForPlotting.csv'),
           sep = ';',
           dec = ',', stringsAsFactors = TRUE)
FirstGroupHDSListDataFrames <- 
  FirstGroupHDSListDataFrames[FirstGroupHDSListDataFrames$conditionGroupPlotting != 'not included',]
# should I take an average per model? 

## Plot first group

plot1 <- ggplot( FirstGroupHDSListDataFrames, 
               aes( x=conditionGroupPlotting , y=score , color=model)) +
  geom_jitter(alpha = 0.5) + scale_color_colorblind() +
  theme( text = element_text(size=11)) +
  labs(x = "grouped condition" , y = "HDS") +
  stat_summary(fun = mean,geom = 'crossbar') 

plot1ColorGroup <- ggplot( FirstGroupHDSListDataFrames, 
                 aes( x=conditionGroupPlotting , 
                      y=score, 
                      color = conditionGroupPlotting)) +
  geom_jitter() + scale_color_colorblind() +
  theme( text = element_text(size=11)) +
  labs(x = "grouped condition" , y = "HDS") +
  stat_summary(fun = mean,geom = 'crossbar') 

plot1Violin <- ggplot( FirstGroupHDSListDataFrames, 
                           aes( x=conditionGroupPlotting , 
                                y=score, 
                                color = conditionGroupPlotting)) +
  geom_violin() + scale_color_colorblind() +
  theme( text = element_text(size=11)) +
  labs(x = "grouped condition" , y = "HDS") +
  stat_summary(fun = mean,geom = 'crossbar') 

## which one to keep? 

# plot second group

SecondGroupHDS <- 
  read.csv(file = paste0(outputPath,
                         'LeaveOneOutCrossValidationSecondGroupForPlotting.csv'),
           sep = ';',
           dec = '.', stringsAsFactors = TRUE)
SecondGroupHDS <- 
  SecondGroupHDS[SecondGroupHDS$conditionGroupPlotting != 'not included',]
## considering leaving fructose out [SecondGroupHDS$GEOid != "GSE123894",]
plot2 <- ggplot(SecondGroupHDS, 
                 aes( x=conditionGroupPlotting , 
                      y=score , 
                      color=model)) +
  geom_jitter(alpha = 0.5) + 
  scale_color_colorblind() +
  theme( text = element_text(size=11)) +
  labs(x = "grouped condition" , y = "HDS") +
  stat_summary(fun = mean,geom = 'crossbar') 

plot2ColorGroup <- ggplot(SecondGroupHDS, 
                           aes( x=conditionGroupPlotting , 
                                y=score, 
                                color = conditionGroupPlotting)) +
  geom_jitter() + scale_color_colorblind() +
  theme( text = element_text(size=11)) +
  labs(x = "grouped condition" , y = "HDS") +
  stat_summary(fun = mean,geom = 'crossbar') 

panelGroups1and2 <- plot1 | plot2

ggsave('GroupedModelsCrossValidation.png', path = outputPath)

# plot acetominophen on it's own because it's a time series: GSE111828 

p3 <- ggplot(HDS_cv.lomo.ListDataFrames[[7]], aes(x = condition,
                                            y = score,
                                            color = condition)) + 
  geom_jitter(alpha = 0.5 ) + scale_color_colorblind() +
  theme( text = element_text(size=11)) +
  labs(x = "time post acetominophen application" , y = "HDS") +
  stat_summary(fun = mean, geom = 'crossbar') 

colorsAcetaminophen <-c('0hr' = 'black', '12hr' = '#E69F00', '24hr' = '#E69F00', 
              '36hr' = '#E69F00', '48hr' = '#E69F00', '72hr' = '#E69F00')

p3TwoColor <- ggplot(HDS_cv.lomo.ListDataFrames[[7]], aes(x = condition,
                                                         y = score,
                                                         color = condition)) + 
  geom_jitter(alpha = 0.5 ) + 
  theme( text = element_text(size=11)) +
  labs(x = "time post acetominophen application" , y = "HDS") +
  stat_summary(fun = mean, geom = 'crossbar') +
  scale_color_manual(values = colorsAcetaminophen)


panel <- (plot1 | plot2) / p3 + plot_annotation(tag_levels = 'A')


ggsave('GroupedModelsANDAcetaminophenCrossValidation.png', path = outputPath)
p3
ggsave('AcetaminophenCrossValidation.png', path = outputPath)

## now aging data set 

levelOrder <- c('1', '3', '6', '9', '12', '15', '18', '21', '24', '27')

HDS_cv.lomo.ListDataFrames[[13]]$condition <- 
  factor(HDS_cv.lomo.ListDataFrames[[13]]$condition, levels = levelOrder)

colorsAge <-c('1' = 'black', '3' = '#E69F00', '6' = '#E69F00', 
              '9' = '#E69F00', '12' = '#E69F00', '15' = '#E69F00',
              '18' = '#E69F00', '21' = '#E69F00', '24' = '#E69F00',
              '27' = '#E69F00')

p4 <- ggplot(HDS_cv.lomo.ListDataFrames[[13]], aes(x = condition,
                                                  y = score,
                                                  color = condition)) + 
  geom_jitter(alpha = 0.5 ) + #scale_color_colorblind() +
  theme( text = element_text(size=11)) +
  labs(x = "age in months" , y = "HDS") +
  stat_summary(fun = mean, geom = 'crossbar') +
  scale_color_manual(values = colorsAge) 


temp <- HDS_cv.lomo.ListDataFrames[[13]]
temp$condition <- as.numeric(as.character(temp$condition))

p4RegressionLineCor <- ggplot(temp, aes(x = condition,
                                        y = score)) + 
  geom_point(alpha = 0.5 ) + 
  theme( text = element_text(size=11)) +
  labs(x = "age in months" , y = "HDS") + 
  stat_smooth(method = 'lm', linetype='dashed', color= '#E69F00') + 
  stat_cor(method = 'spearman', label.x = 3, label.y = -0.1 ) + 
  scale_x_discrete(limits=c(1, seq(3,27,3)))


p4
ggsave('AgeCrossValidation.png', path = outputPath)
panelVersion3 <-(plot1 | plot2) / p3TwoColor / p4 + 
  plot_annotation(tag_levels = 'A')
panelVersion3
ggsave('AllResutsCrossValidationManyColors.png', path = outputPath)

p4RegressionLineCor
ggsave('AgeCrossValidationWithLM.png', path = outputPath)

panelVersion4 <-(plot1 | plot2) / p3TwoColor / p4RegressionLineCor + 
  plot_annotation(tag_levels = 'A')
panelVersion4
ggsave('AllResutsCrossValidationManyColorsRegress.png', path = outputPath)


# OLD APPROACH
# create list of plots: there's 13 plots 
## relevel datt$condition <- relevel( datt$condition , ref = "StSt")

HDS_cv.lomo.ListBoxplots <- lapply( seq(HDS_cv.lomo.ListDataFrames) , function(ii){
  ppg <- ggplot( HDS_cv.lomo.ListDataFrames[[ii]] , aes( x=condition , y=score , color=condition)) +
    geom_jitter(shape = 2, size = 3) + theme_bw() + theme( text = element_text(size=18)) +
    labs(x = "condition" , y = "HDS") +  stat_compare_means() + 
    stat_summary(fun = mean,geom = 'crossbar')
  
  return( ppg )
}
)

HDS_Boxplots <- lapply( seq(HDS_cv.lomo.ListDataFrames) , function(ii){
  ppg <- ggplot( HDS_cv.lomo.ListDataFrames[[ii]] , 
                 aes( x=condition , y=score , color=condition)) +
    geom_jitter(shape = 2, size = 3, stroke = 1) + 
    scale_color_colorblind() +
    theme( text = element_text(size=11)) +
    labs(x = "condition" , y = "HDS") 
  
  return( ppg )
}

 
)  
names(HDS_Boxplots) <- names(HDS_cv.lomo)

panel <-  HDS_Boxplots[["EtOH"]] | HDS_Boxplots[["STAM"]] | HDS_Boxplots[["DEN"]]
  
