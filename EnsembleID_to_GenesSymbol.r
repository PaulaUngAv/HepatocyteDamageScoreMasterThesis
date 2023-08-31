library(biomaRt)

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
searchDatasets(mart = ensembl, pattern = "mmusculus")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

# CASE 1: Ensemble IDs instead of Gene Symbols 

count_table <- read.delim('GSE119441/GSE119441_genecounts.txt')
tail(count_table,20)
#get rid of rows with p-value = NA 
#count_table <- count_table[is.na(count_table$pvalue) == FALSE ,]
genes <- count_table$id
count_table$Gene.Symbol <- NA
# for transcript ensembl IDs do this: 
genes <- sub( "\\..*","",genes)
count_table$id <- sub('\\..*','', count_table$id ) 
#searchFilters(mart = ensembl, pattern = "symbol")
#searchFilters(mart = ensembl, pattern = "ensemble*")
G_list <- getBM( filter = 'ensembl_gene_id', 
                 attributes= c("ensembl_gene_id","mgi_symbol"), 
                 values = genes,mart= ensembl)

# regard that, when ensembl ID has not matching gene symbol in G_list, 
# the following command will assign NA as the name of the gene 
count_table$Gene.Symbol <- G_list$mgi_symbol[match( count_table$id,
                                                    G_list$ensembl_gene_id 
)]
#this shows us how many genes could not get assigned a gene symbol
length(count_table$Gene.Symbol[is.na(count_table$Gene.Symbol)] )

# remove them Gene.Symbols == NAs 
count_table <- count_table[is.na(count_table$Gene.Symbol) == FALSE ,]

#does this match: YES -> bug fixed 
length( unique( G_list$mgi_symbol))
length( unique( count_table$Gene.Symbol))

# this does not make sense 
count_table.ag.p <- aggregate(cbind(pvalue, log2FoldChange) ~ Gene.Symbol, 
                              data = count_table, FUN = mean )

#check if final length makes sense
length( unique( count_table.ag.p$Gene.Symbol))

# remove first row if empty
count_table <- count_table.ag.p[-1,]


## save 
write.csv(count_table_test, 
          'GSE119441/GSE119441_count_table.csv', 
          row.names = F)

# CASE 2: transcript ensembl IDs and gene symbol available
## this means, some genes might come up twice and you dont need to convert 
# to gene symbol but make sure you aggregate different transcripts 

# GSE135050 _hfcfdiet_deg
# here, gene symbols included so just aggregate the repeated genes 
count_table <- read.delim('Desktop/Project Module Beyer /GSE135050/GSE135050 _hfcfdiet_deg_ensemlids.csv', 
                          sep = ',')
#get rid of rows with p-value = NA 
count_table <- count_table[is.na(count_table$pvalue) == FALSE ,]
#p-adjusted has many NAs and messes up table
count_table.ag.p <- aggregate(cbind(pvalue, log2FoldChange) ~ Gene.Symbol, 
                              data = count_table, FUN = mean )

#compare to see if it makes sense
length( unique( count_table$Gene.Symbol))
length( unique( count_table.ag.p$Gene.Symbol))

#if first row empty
#count_table <- count_table.ag.p[,-1]
count_table <- count_table.ag.p

write.csv(count_table, 'Desktop/Project Module Beyer /GSE135050/GSE135050_hfcfdiet_deg.csv',row.names = F)
