# In Progress, do not use

doTTest = function() {
    read.table("Desktop/R Project Desktop/full_BRCA_miR_data.txt",
               stringsAsFactors = F, row.names = 1, header = T) -> miRNA
    
    for (i in nrows(miRNA)){
        curRow = as.numeric(miRNA[i,])
        for (j in nrows(miRNA)-i)
    }
    
}