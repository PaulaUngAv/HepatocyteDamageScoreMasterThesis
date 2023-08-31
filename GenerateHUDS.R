## Hepatocyte Damage Score 
library(ggplot2)
library(ggpubr)
library(Seurat)
#### READ ME 
# script to create latest version of HUDS

# load 1) function that generates HUDS
source(file='~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R')

path <- '~/Desktop/HepatocyteDamageScore/'
# directory containing DEG tables 
setwd('~/Desktop/HepatocyteDamageScore/All_DEG_Tables/')
# total of 18 DEG tables, no KO studies
studies <- c('GSE99010_CCl4_deg.csv',
             'GSE99010_WesternDiet_deg.csv',
             'GSE97234_AH_deg.csv',
             'GSE97234_ASH_deg.csv',
             'GSE83240_DEN_deg.csv',
             'GSE153580_STZ_deg.csv',
             'GSE148849_fastfooddiet_deg.csv',
             'GSE138419_AMLN_deg.csv',
             'GSE137449_CDAHFD_diet_deg.csv',
             'GSE135050_hfcfdiet_deg.csv',
             'GSE132040_young(6_9_12)vs_old(18_21_24)_deg.csv',
             'GSE123894_Fructose_only_DBA2Jstrain_deg.csv',
             'GSE119953_HFD_deg.csv',
             'GSE119953_DDC_deg.csv',
             'GSE119441_HFD_deg.csv',
             'GSE119441_PFOA_deg.csv',
             'GSE114261_STAM_deg.csv',
             'GSE111828_Acetaminophen_deg.csv',
             'GSE97024_Oxr_BrainKO_deg.csv',
             'GSE61117_Trim24_KO_deg.csv',
             'GSE124694_Arid1a_LiverKO_deg.csv',
             'GSE109431_LHR1_KO_deg.csv',
             'GSE106369_MOFKO_deg.csv'
)

DElist <- lapply( seq( studies ), function(ii){
  read.csv(studies[[ii]])
})

setwd(path)

# Load list of genes expressed in hepatocytes
pathGeneFilter <- '~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/'
# I have three geneFilter lists, one from back then:
# geneFilter <- read.csv(file =  
#                          paste0(pathGeneFilter,
#                                 'p_nucseq_genes_filtered_genelist_minus5threshold_log2.csv'))
# # 2) "my method" genes expressed in 25% or more nuclei sn data set of Liver Cell Atlas
 geneFilter <- read.csv(file ='~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/p_nucseq_genes_filtered_genelist_minus5threshold_log2TESTtoCOMPARE.csv')
# 3) Tim's method: genes expressed in > 0.5% nuclei sn Liver Cell Atlas data
#geneFilter <- read.csv(file =
#                         paste0(pathGeneFilter,
#                                'GeneListTimsMethod.csv'))

geneFilter <- geneFilter$Genes

# Create HUDS and store as csv file
final_table <- damage_signature.func(DElist, geneFilter = geneFilter)

# HUDS from 1) original method
# write.csv(final_table, 'Universal_Hepatocyte_Damage_Signature_9June22.csv')

final_table <- final_table[order(final_table$mean_rank, decreasing = FALSE),] 
rownames(final_table) <- 1:length(final_table$gene_symbol)

# 2) HUDS made with Tim's method gene filter
write.csv(final_table, paste0('~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS_GeneFilterp_nucseq_genes_filtered_genelist_minus5threshold_log2TESTtoCOMPARE.csv'))
# test to compare: 
geneSetNow <- final_table$gene_symbol
genesSetThen <- read.csv(file = '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/Universal_Hepatocyte_Damage_Signature_9June22.csv')
genesSetThen <- genesSetThen$gene_symbol

length(intersect(head(geneSetNow, 100), head(genesSetThen, 100)))
length(intersect(geneSetNow, genesSetThen))
length(genesSetThen)
length(geneSetNow)
