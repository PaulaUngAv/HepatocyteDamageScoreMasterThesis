# Generate Hepatocyte Universal Damage Signature #

# load functions written by us
# generate damage signature
source("~/Repositories/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
# calculate damage score
source("~/Repositories/cell-damage-score/AUCell_script.r")

# load DEG tables
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

DElist <- lapply( seq( studies ), function(ii){
  read.csv(paste0(path,"AllBulkStudiesDEGAnalysis/",studies[[ii]]))
})

names(DElist) <- basename(studies)

# laod list of genes to be considered as possibly expressed in hepatocytes
geneFilter <- 
  read.csv(file = 
             '~/Desktop/HepatocyteDamageScore/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv',
           row.names = 1)
geneFilter <- geneFilter$Genes

# create HUDS

HUDS <- damage_signature.func(DElist = DElist, 
                      geneFilter = geneFilter)

rownames(HUDS)<- 1:length(HUDS$gene_symbol)

write.csv(HUDS,'~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')

