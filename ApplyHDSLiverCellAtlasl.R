library(Seurat)
library(AUCell)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)

#### load functions written by us
# generate damage signature
{
  source("~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
  # calculate damage score
  source("~/Repositories/cell-damage-score/AUCell_script.r")
}

outputPath <- '~/Desktop/HepatocyteDamageScore/Results/SingleNucleiCell/'


# Hepatocyte Universal Damage Signature
HUDS <- read.csv(file = 
                   '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')

# get counts 
data <- readRDS(
  '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv')
annotationPath <- 
  '~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data/ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv'
annotation <- read.csv(annotationPath, row.names = 1) 
sampleInfo <- read.csv('~/Desktop/HepatocyteDamageScore/scRNAseq/Liver Cell Atlas Data//GSE192740_sampleInfo_modified.csv', sep = ";")
HDS <- DS_calc.func(data, DSignature = HUDS,
                    useOrder = 'mean_rank',
                    ntop = 42)


HDSdataframe <- data.frame('Cell_ID' = names(HDS),
                     'DS' = unname(HDS))

HDSdataframe$diet[match(annotation$cell_id, HDSdataframe$Cell_ID)] <- annotation$condition
HDSdataframe$sample[match(annotation$cell_id, HDSdataframe$Cell_ID)] <- annotation$Sample

## diet here is not the acutal DIET of the mice but the data set they belong to!! 
HDSdataframe$diet <- factor(HDSdataframe$diet,
                            ordered = TRUE,
                            levels = c("StSt", "Nafld"))

HDSdataframe$Annotation <- NA

HDSdataframe[HDSdataframe$sample == "CS197" | 
               HDSdataframe$sample ==  "CS196" |
               HDSdataframe$sample == "CS193" |
               HDSdataframe$sample == "CS192" |
               HDSdataframe$sample == "CS191" |
               HDSdataframe$sample == "CS185"  |
               HDSdataframe$sample == "CS183" , 'Annotation'] <- 'Western Diet'

HDSdataframe$Annotation[is.na(HDSdataframe$Annotation)] <- 'Standard Diet'

HDSdataframe$Annotation <- factor(HDSdataframe$Annotation,
                            ordered = TRUE,
                            levels = c("Standard Diet", "Western Diet"))

# add age

HDSdataframe$age <- NA
HDSdataframe[HDSdataframe$sample == "CS197" | 
               HDSdataframe$sample ==  "CS196" |
               HDSdataframe$sample == "CS195" |
               HDSdataframe$sample == "CS194" , 'age'] <- '36 weeks'

HDSdataframe[HDSdataframe$sample == "CS193" | 
               HDSdataframe$sample ==  "CS192" |
               HDSdataframe$sample == "CS191" |
               HDSdataframe$sample == "CS185" |
               HDSdataframe$sample == "CS183" |
               HDSdataframe$sample == "CS190" |
               HDSdataframe$sample == "CS189" |
               HDSdataframe$sample == "CS188" |
               HDSdataframe$sample == "CS184" |
               HDSdataframe$sample == "CS182", 'age'] <- '24 weeks'

HDSdataframe$age <- factor(HDSdataframe$age,
                                  ordered = TRUE,
                                  levels = c('24 weeks', '36 weeks'))

saveRDS(HDSdataframe, file = paste0(outputPath, 'HDSLiverCellAtlasAllHep.rds'))

## Plot results of Nafld data set 
plotNafldSampleFacetAges <- ggplot( HDSdataframe[HDSdataframe$diet == "Nafld",], 
                        aes( x = Annotation, 
                             y = DS,
                             color = Annotation
                        )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind() +
  theme( text = element_text(size=12),
         plot.title = element_text(size=12),
         legend.position = 'right',
         axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5
                                    )) +
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Diet"))+
  labs(y = "HDS", x = '') + 
  coord_cartesian(ylim = c(-0.5, 0.4)) + 
  scale_y_continuous(breaks=c(-0.4, -0.2, 0, 0.2, 0.4)) +
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") +
  facet_wrap(~age)

## plot distribution of data set healthy mouse (ABU samples)
plotBatchSt <- ggplot( HDSdataframe[HDSdataframe$diet == "StSt",], 
                        aes( x = diet, 
                             y = DS,
                             color = diet
                        )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind() +
  theme( text = element_text(size=12),
         plot.title = element_text(size=12),
         legend.position = 'none') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS", x = '') + 
  coord_cartesian(ylim = c(-0.5, 0.4)) + 
    scale_y_continuous(breaks=c(-0.4, -0.2, 0, 0.2, 0.4)) +
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") 

# not included in thesis 
# plotAll <- ggplot( HDSdataframe, 
#                    aes( x = Annotation, 
#                         y = DS,
#                         color = Annotation
#                    )) + 
#   geom_violin(trim = FALSE) + 
#   scale_color_colorblind() +
#   theme( text = element_text(size=12),
#          plot.title = element_text(size=12),
#          legend.position = 'none',
#          axis.text.x = element_text(angle = 45, 
#                                     vjust = 0.5, hjust=0.5
#          )) +
#   guides(colour = guide_legend(override.aes = aes(label = ""),
#                                title="Diet"))+
#   labs(y = "HDS", x = '') + 
#   coord_cartesian(ylim = c(-0.5, 0.4)) + 
#   scale_y_continuous(breaks=c(-0.4, -0.2, 0, 0.2, 0.4)) +
#   stat_summary(fun.data = "mean_sdl",
#                fun.args = list(mult = 1),
#                geom = "pointrange")


# panel 

panelMA <- (plotBatchSt | plotNafldSampleFacetAges) + 
  plot_annotation(tag_levels = 'A')

ggsave('LiverCellAtlasAllSingleNucleiHDS.png', path = outputPath)

## summary stats
### NAFLD data set
## 24 weeks 
# st diet
summary(HDSdataframe[HDSdataframe$diet == "Nafld" & 
                       HDSdataframe$age == '24 weeks' & 
                       HDSdataframe$Annotation != 'Western Diet',])
# wd diet
summary(HDSdataframe[HDSdataframe$diet == "Nafld" & 
                       HDSdataframe$age == '24 weeks' & 
                       HDSdataframe$Annotation == 'Western Diet',])

## 35 weeks 
# st diet
summary(HDSdataframe[HDSdataframe$diet == "Nafld" & 
                       HDSdataframe$age == '24 weeks' & 
                       HDSdataframe$Annotation != 'Western Diet',])
# wd diet
summary(HDSdataframe[HDSdataframe$diet == "Nafld" & 
                       HDSdataframe$age == '24 weeks' & 
                       HDSdataframe$Annotation == 'Western Diet',])


# StSt data set
summary(HDSdataframe[HDSdataframe$diet == "Nafld",])

