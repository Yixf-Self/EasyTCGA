
#' Reshape sample-level log2 miRSeq expression values of the sample type TP. Primarily this service is intended to reshape all sample-level log2 miRSeq expression values of one cohort.
#' See \code{\link{dn_miRSeq_cohort}}.
#' @param data data.frame of log2 miRSeq expression values. See \code{\link{dn_miRSeq}} for downloading these values.
#' @return nxp-matrix, n = number of patients, p = number of mirs, (matrix(i,j))_{i,j} = ((data$tcga_participant_barcode==barcode[i], data$mir==mir[j]))_{i,j}
#' @export
reshape.miRSeq = function(data, sample_type = "TP"){

  cat("Info: the sample types of the data are:", unique(data$sample_type) )
  data = subset(data, data$sample_type=="TP")
  idx = which(colnames(data)=='expression_log2')
  barcode = unique(data[,"tcga_participant_barcode"])
  mir = unique(data[,"mir"])
  n = length(barcode)
  p = length(mir)
  matrix = matrix(0, nrow = n, ncol = p)
  colnames(matrix) = mir
  rownames(matrix) = barcode

  # fill in column by column:
  for (j in 1:p) {
    matrix[,j] =  as.numeric(data[(data$mir==mir[j]), idx])
  }

  return(matrix)

}




#' Reshape sample-level log2 mRNASeq expression values of the sample type TP. Primarily this service is intended to reshape all sample-level log2 mRNASeq expression values of one cohort.
#' See \code{\link{dn_mRNASeq_cohort}}.
#' @param data data.frame of log2 mRNASeq expression values. See \code{\link{dn_mRNASeq}} for downloading these values.
#' @return nxp-matrix, n = number of patients, p = number of genes, (matrix(i,j))_{i,j} = ((data$tcga_participant_barcode==barcode[i], data$gene==gene[j]))_{i,j}
#' @export
reshape.mRNASeq = function(data, sample_type = "TP"){

  cat("Info: the sample types of the data are:", unique(data$sample_type))
  data = subset(data, data$sample_type=="TP")
  idx = which(colnames(data)=='expression_log2')
  barcode = unique(data[,"tcga_participant_barcode"])
  gene = unique(data[,"gene"])
  n = length(barcode)
  p = length(gene)
  matrix = matrix(0, nrow = n, ncol = p)
  colnames(matrix) = gene
  rownames(matrix) = barcode

   # fill in column by column:
  for (j in 1:p) {
    matrix[,j] =  as.numeric(data[(data$mir==mir[j]), idx])
  }

  return(matrix)

}

#'
#'
#'
#' #' @examples
#' #' cohort = "ESCA"
#' #' page_size = 2000
#' #' esca.miRSeq = dn_miRSeq_cohort(cohort, page.Size)
#' #' esca.miRSeq_reshaped = reshape.miRSeq (esca.miRSeq, sample_type = "TP")
#' #' esca.mRNASeq = dn_mRNASeq_cohort(cohort, page.Size)
#' #' esca.mRNASeq_reshaped = reshape.mRNASeq (esca.mRNASeq, sample_type = "TP")
#'
#'
#'
