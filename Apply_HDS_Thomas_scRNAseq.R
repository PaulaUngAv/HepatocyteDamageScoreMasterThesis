# Applying HDS to Thomas scRNAseq data
library(Seurat)
library(AUCell)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gridExtra)

data <- readRDS('~/Desktop/HepatocyteDamageScore/scRNAseq/Thomas_Wunderlich_snRNAseq/nik_samples_annotated.rds')
data <- subset(data, subset= class_annotation == "Hepatocyte")
metadata_hepatocytes <- data@meta.data

data_nash <- subset(data, subset = Diet == 'Nash')
data_cdaa <- subset(data, subset = Diet == 'Cdaa')
data_hcc <- subset(data, subset = Diet == 'Hcc')

remove(data)

#### load functions
# generate damage signature
source("~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
# calculate damage score
source("~/Repositories/cell-damage-score/AUCell_script.r")

#DSignature
UniHepDamageSignature <- 
  read.csv(file = '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')

DS_nash <- DS_calc.func( exprMatrices = GetAssayData(object = data_nash,
                                                     slot = 'counts'), 
                         ceilThrsh = 0.05 , 
                         DSignature = UniHepDamageSignature , 
                         ntop = 42, useOrder ="mean_rank" )

DS_cdaa <- DS_calc.func( exprMatrices = GetAssayData(object = data_cdaa,
                                                     slot = 'counts'), 
                         ceilThrsh = 0.05 , 
                         DSignature= UniHepDamageSignature , 
                         ntop = 42, useOrder ="mean_rank" )

DS_hcc <- DS_calc.func( exprMatrices = GetAssayData(object = data_hcc,
                                                    slot = 'counts'), 
                         ceilThrsh = 0.05 , 
                         DSignature= UniHepDamageSignature , 
                         ntop = 42, useOrder ="mean_rank" )


#Format data for plotting

DS_all <- data.frame('Cell_ID' = names(DS_nash),
                     'DS' = unname(DS_nash),
                     'Diet' = rep('NASH', length(names(DS_nash))))

DS_all <- rbind(DS_all, data.frame('Cell_ID' = names(DS_cdaa),
                                   'DS' = unname(DS_cdaa),
                                   'Diet' = rep("CDAA", length(names(DS_cdaa)))))

DS_all <- rbind(DS_all, data.frame('Cell_ID' = names(DS_hcc),
                                   'DS' = unname(DS_hcc),
                                   'Diet' = rep("HCC", length(names(DS_hcc)))))

DS_all$Genotype <- unlist(lapply(DS_all$Cell_ID, function(ii){
  return(metadata_hepatocytes[match(ii, metadata_hepatocytes$Cell_ID), 'Genetic']) 
}))

DS_all$Genotype <- factor(DS_all$Genotype,
                      ordered = TRUE,
                      levels = c("NIK-FL", "NIK-KO"))

#ploting 

plotNASH <- ggplot( DS_all[DS_all$Diet == 'NASH',], 
        aes( x = Genotype  , y = DS,
             color = Genotype)) + 
  geom_violin(trim = FALSE) + #theme_bw() + 
  scale_color_colorblind() +
  ggtitle("NASH") +
  theme( text = element_text(size=12),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(size=12),
         legend.position = 'none') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS") +
  coord_cartesian(ylim = c(-0.5, 0.5)) + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange"
  ) 

plotCDAA <- ggplot( DS_all[DS_all$Diet == 'CDAA',], 
                    aes( x = Genotype  , y = DS,
                         color = Genotype)) + 
  geom_violin(trim = FALSE) + #theme_bw() + 
  scale_color_colorblind() +
  ggtitle("CDAA") +
  theme( text = element_text(size=12),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(size=12),
         legend.position = 'none') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS") +
  coord_cartesian(ylim = c(-0.5, 0.5)) + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") 

plotDEN <- ggplot( DS_all[DS_all$Diet == 'HCC',], 
                    aes( x = Genotype  , y = DS,
                         color = Genotype)) + 
  geom_violin(trim = FALSE) + #theme_bw() + 
  scale_color_colorblind() +
  ggtitle("DEN") +
  theme( text = element_text(size=12),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(size=12),
         legend.position = 'right') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS") + 
  coord_cartesian(ylim = c(-0.5, 0.5)) + 
  scale_y_continuous(breaks=c(-0.5,-0.25, 0, 0.25, 0.5)) +
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") 

panelWunderlich <- (plotNASH | plotCDAA | plotDEN) + 
  plot_annotation(tag_levels = 'A')

ggsave(filename = 'HDSdistributionNIKallHeps.png', 
       path = '~/Desktop/HepatocyteDamageScore/Results/SingleNucleiCell'
       )

## summary statistics

summary(DS_all)

# NASH
summary(DS_all[DS_all$Diet == 'NASH' &
                 DS_all$Genotype == 'NIK-FL', ])

summary(DS_all[DS_all$Diet == 'NASH' &
                 DS_all$Genotype == 'NIK-KO', ])

# CDAA 
summary(DS_all[DS_all$Diet == 'CDAA' &
                 DS_all$Genotype == 'NIK-FL', ])
summary(DS_all[DS_all$Diet == 'CDAA' &
                 DS_all$Genotype == 'NIK-KO', ])

# DEN
summary(DS_all[DS_all$Diet == 'HCC' &
                 DS_all$Genotype == 'NIK-FL', ])
summary(DS_all[DS_all$Diet == 'HCC' &
                 DS_all$Genotype == 'NIK-KO', ])

