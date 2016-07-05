

reshape.miRSeq = function(data){
  
  idx = which(colnames(data)=='expression_log2') # index of column containing log2 miRSeq expression values
  barcode = unique(data[,"tcga_participant_barcode"]) 
  mir = unique(data[,"mir"]) 
  n = length(barcode)
  p = length(mir)
  matrix = matrix(, nrow = n, ncol = p)
  colnames(matrix) = mir
  rownames(matrix) = barcode
  cat(dim(matrix), "")
  # fill in column by column:
  for (j in 1:p) {
    cat(j,"")
    for (i in 1:n) {
      # choose i and j so that matrix(i,j) = expression_log2(barcode(i), mir(j))
      # values of the column j
      matrix[i,j] =  data[(data$mir==mir[j]) & (data$tcga_participant_barcode==barcode[i]),idx] 
     cat(matrix[i,j], "")
      cat("", i, "")
    }
  # cat(matrix[80,1])
}

return(matrix)

}


reshape.mRNASeq = function(data){
  
  idx = which(colnames(data)=='expression_log2') # index of column containing log2 mRNASeq expression values
  barcode = unique(data[,"tcga_participant_barcode"]) 
  gene = unique(data[,"gene"]) 
  n = length(barcode)
  p = length(gene)
  matrix = matrix(, nrow = n, ncol = p)
  colnames(matrix) = gene
  rownames(matrix) = barcode
  cat(dim(matrix), "")
  # fill in column by column:
  for (j in 1:p) {
    cat(j,"")
    for (i in 1:n) {
      # choose i and j so that matrix(i,j) = expression_log2(barcode(i), mir(j))
      # values of the column j
      matrix[i,j] =  data[(data$gene==gene[j]) & (data$tcga_participant_barcode==barcode[i]),idx] 
      cat(matrix[i,j], "")
      cat("", i, "")
    }
    # cat(matrix[80,1])
  }
  
  return(matrix)
  
}


reshape.mutaton = function(data){
  
  idx = which(colnames(data)=='Protein_Change') # index of column containing protein change values
  barcode = unique(data[,"tcga_participant_barcode"]) 
  Hugo_Symbol = unique(data[,"Hugo_Symbol"]) 
  n = length(barcode)
  p = length(Hugo_Symbol)
  matrix = matrix(, nrow = n, ncol = p)
  colnames(matrix) = Hugo_Symbol
  rownames(matrix) = barcode
  cat(dim(matrix), "")
  # fill in column by column:
  for (j in 1:p) {
    cat(j,"")
    for (i in 1:n) {
      # choose i and j so that matrix(i,j) = expression_log2(barcode(i), mir(j))
      # values of the column j
      matrix[i,j] =  data[(data$Hugo_Symbol==Hugo_Symbol[j]) & (data$tcga_participant_barcode==barcode[i]),idx] 
      cat(matrix[i,j], "")
      cat("", i, "")
    }
    # cat(matrix[80,1])
  }
  
  return(matrix)
  
}