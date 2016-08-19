# "Imports" necessary files
source("MergeTables.R")


# Change these values
mRNA = read.table("Trimmed_BRCA_mRNA.txt", header=T, row.names=1, stringsAsFactors=F)
miRNA = read.table("GroupB_BRCA_miR.txt", header=T, row.names=1, stringsAsFactors=F)
fullTable = mergeTablesByPatientID(mRNA, miRNA)
partialTable = as.matrix(splitmiRNATable(as.matrix(fullTable)))




corList = c()

# iter1 is the index of the first part of the pair
# iter2 is the index of the row being checked against
for (iter1 in 1:(nrow(partialTable)-1)) {       # iterates over each row in the table, except for the last (right before the last index is reached, all pairs will have been checked, and errors may ensue)
    for (iter2 in (iter1+1):(nrow(partialTable))) { # iter1+1 prevents comparing a row against itself (which would be bad)
        # Adds the correlation of the current pair to a vector
        # use="complete" prevents failing on NA
        
        corList = c(corList,
                    cor(partialTable[iter1,], partialTable[iter2,]))

        # Name the correlation value just added.
        # Without this, it would be impossible to match sorted correlations to the miRNA pairs.
        names(corList)[length(corList)] =   # Gets the last element of the vector; the last element added
            paste(row.names(partialTable)[iter1], row.names(partialTable)[iter2], sep=" vs. ")
       
    }
}
corList = rev(sort(abs(corList))) # Lists from most highly correlated to lowest
highCor=corList[1:100] # Top 100 most highly correlated in vector called "highCor"
highCor[1:10]
