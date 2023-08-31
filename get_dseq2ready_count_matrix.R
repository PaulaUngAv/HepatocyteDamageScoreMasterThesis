# Create Count Matrix 
#set working directory to the one that contains all reads 
setwd(dir = 'GSE83240/')
run_list <- read.delim('SRR_Acc_List.txt',header = F)
run_list <- c(run_list[,1])
run_list <- run_list[order(run_list)]
add <- 'ReadsPerGene.out.tab'
filename <- paste0(run_list[1],add)

# initiate count table with first run
temp <- read.delim(filename)
temp <- temp[-1:-3,]
count_table <- data.frame('gene_symbol' = temp[,1], 'run1' = temp[,2])
names(count_table)[names(count_table) == 'run1'] <- run_list[1]

# add the other runs
i = 2
for (i in 2:(length(run_list))) {
  run <- run_list[i]
  filename <- paste0(run,add)
  print(filename)
  temp <- read.delim(file = filename)
  temp <- temp[-1:-3,]
  
  if (identical(count_table[,1],temp[,1]) == TRUE){
    count_table$temp <- temp[,2]
    names(count_table)[names(count_table) == 'temp'] <- run
    print('yay!')
  } else { print('error: different genes')}
  
}


write.table(count_table, file = 'GSE83240_count_table', sep = ' ', row.names = F,
            col.names = T)
remove(count_table)
remove(temp)
