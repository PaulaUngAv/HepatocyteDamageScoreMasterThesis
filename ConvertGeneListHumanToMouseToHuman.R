##########################################################
### Convert mouse gene list to human gene list and vs. ###
##########################################################
# Tim's original function modified by Paula 

### function to convert Human to mouse genes
genes <- read.csv('~/Desktop/HepatocyteDamageScore/Universal Damage Signature/Universal_Hepatocyte_Damage_Signature_9June22.csv')
genes <- genes$gene_symbol

fun_homoTO.FROMmouse <- function(genes, TO = TRUE ){
  # if TO parameter is TRUE then Homo -> Mouse, if FALSE then Mouse -> Homo
  
  options(connectionObserver = NULL)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(Orthology.eg.db)
  
  if(isTRUE(TO)){
    egs <- mapIds(org.Hs.eg.db, genes, "ENTREZID","SYMBOL")
    mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
    mapped$MUS <- mapIds(org.Mm.eg.db, 
                         as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
    mapped$MUS.ens <- mapIds(org.Mm.eg.db, 
                             as.character(mapped$Mus.musculus), "ENSEMBL", "ENTREZID")
  } else {
    egs <- mapIds(org.Mm.eg.db, genes, "ENTREZID","SYMBOL")
    mapped <- select( Orthology.eg.db, egs, "Homo.sapiens", "Mus.musculus")
    mapped$HOMO <- mapIds(org.Hs.eg.db, 
                          as.character(mapped$Homo.sapiens), "SYMBOL", "ENTREZID")
    mapped$HOMO[ is.na(mapped$Homo.sapiens)] <- NA
     # toupper(rownames(mapped))[ is.na( mapped$Homo.sapiens ) ]
    mapped$HOMO[rownames(mapped)=="Cd59a"]<- "CD59"
    mapped$HOMO.ens <- mapIds(org.Hs.eg.db, 
                              as.character(mapped$HOMO), "ENSEMBL", "SYMBOL")
    
  }
  return(mapped)
}

mapped_MusToHom <- fun_homoTO.FROMmouse(genes, FALSE)
# !!!!! WARNING !!!!!!!!
# beware that genes in column "HumanEnsemble" marked as NA, are not identified
# correctly in "HumanGeneSymbol" column
TopHDSGenesConverted <- 
  data.frame("MouseGeneSymbol" = head(genes,42),
             "HumanEnsemble" = head(mapped_MusToHom$HOMO.ens, 42),
             "HumanGeneSymbol" = head(unname(unlist(mapped_MusToHom$HOMO)), 42))
sum(which(is.na(mapped_MusToHom$Homo.sapiens)) <= 42)
# ==> 10 from 42 top genes not found 
# This will be identified looked into manually and decided what to do about them

write.csv(TopHDSGenesConverted, 
            file = 
              '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDSforHDSMouseToHumanGenesBEFOREManualCorrection.csv')

## generate human version of HUDS leaving ambiguous genes out (marked with NA)

mapped_MusToHom <- head(mapped_MusToHom[!is.na(mapped_MusToHom$Homo.sapiens),],)

TopHDSGenesConverted <- 
  data.frame("MouseGeneSymbol" = row.names(mapped_MusToHom),
             "MouseEntrezID" = mapped_MusToHom$Mus.musculus,
             "HumanEnsembl" = mapped_MusToHom$HOMO.ens,
             "HumanGeneSymbol" = unname(unlist(mapped_MusToHom$HOMO)))
sum(is.na(TopHDSGenesConverted$MouseEntrezID))

write.csv(TopHDSGenesConverted, 
          file = 
            '~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDSforHDSMouseToHumanGenes.csv')


#### Extra section 

# Section to explore on whole list, which genes not matching  
length(unique(mapped_MusToHom$Mus.musculus))
# 728
length(unique(mapped_MusToHom$HOMO.ens))
# 671 ==> 57 genes in HUDS don't have a matching partner 

# mouse genes not mapped to a Homo s. ensemble ID
genesMissing <- unname(unlist(mapped_MusToHom$HOMO[which(is.na(mapped_MusToHom$HOMO.ens))]))
genesMissing <- unname(unlist(mapped_MusToHom$HOMO[which(is.na(mapped_MusToHom$HOMO.ens))]))
# [1]  59 339 356 550 632
# mouse: what genes are not being identified
which(is.na(mapped_MusToHom$Mus.musculus))
# [1]  59 339 356 550 632

