# Collection of functions used to process data from tables

mergeTablesByPatientID <- function(x, y, IDStart = 9, IDEnd = 12) {
    # Merges tables so that same column name are 
    # x and y should be tables with column names as the sample IDs.
    # IDStart and IDEnd are used to split a string, isolating the patient ID.
    
    colX = colnames(x)      # Get names of columns
    colY = colnames(y)      #   and save to new vectors
    
    # Get patient IDs that are in common between the two tables.
    commonX = colX[substr(colX,IDStart,IDEnd)%in%substr(colY,IDStart,IDEnd)]
    commonY = colY[substr(colY,IDStart,IDEnd)%in%substr(colX,IDStart,IDEnd)]
    
    mergedX = x[,commonX]   # Save all the columns that have
    mergedY = y[,commonY]   #   counterparts in the other table.
    
    colnames(mergedX) = substr(commonX,IDStart,IDEnd) # Change the names of the columns so that
    colnames(mergedY) = substr(commonY,IDStart,IDEnd) #   they consist only of the Patient ID.
    
    merged = rbind(mergedX,mergedY)   # Merge the tables
    
    return(merged)    # Return the final result
}


splitmiRNATable = function(x) {
    # x is a matrix processed by mergeTablesByPatientID()
    # This function returns the table with only the miRNA section
    
    return(x[grep("hsa", rownames(x)),])
}


getCorrelated = function(file1, file2) {
    

    mRNA = read.table(file1, header=T, row.names=1, stringsAsFactors=F)
    miRNA = read.table(file2, header=T, row.names=1, stringsAsFactors=F)
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
    corList = rev(sort(corList))
    return(corList)

}
