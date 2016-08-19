# Collection of functions used to process data from tables


mergeTablesByPatientID <- function(x, y, IDStart = 9, IDEnd = 12) {
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

cleanNAs = function(inTable) {
    output = matrix()
    for (i in nrow(inTable)) {
        if () { # Insert code to check if all are NA
            output = rbind(output, inTable[i,])
        }
    }
}