
#' Reshape sample-level log2 miRSeq expression values of the sample type TP. Primarily this service is intended to reshape all sample-level log2 miRSeq expression values of one cohort.
#' See \code{\link{dn_miRSeq_cohort}}.
#' @param data data.frame of log2 miRSeq expression values. See \code{\link{dn_miRSeq}} for downloading these values.
#' @return nxp-matrix, n = number of patients, p = number of mirs, \cr
#'         (matrix(i,j))_{i,j} =  (expression_log2)_{ij}, where \cr
#'         data$tcga_participant_barcode==barcode[i], data$mir==mir[j], \cr
#'         i=1,..,n, j=1,...,p
#' @export  
#' @examples  
#' mir = miRNA_ID[1:10]
#' cohort = "BLCA"
#' tcga_participant_barcode = "TCGA-ZF-AA53"  # TCGA patient barcode from BLCA
#' page.Size = 250
#' sort_by = "mir"
#' obj = dn_miRSeq(mir, cohort, tcga_participant_barcode, page.Size, sort_by)
#' 
#' cohort = "ESCA"     
#' page_size = 2000
#' esca.miRSeq = dn_miRSeq_cohort(cohort, page.Size2)
#' esca.miRSeq_reshaped = reshape.miRSeq (esca.miRSeq, sample_type = "TP")
reshape.miRSeq = function(data, sample_type = "TP"){

  cat("Info: the sample types of the data:", unique(data$sample_type) )
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




#' Reshape sample-level log2 mRNASeq expression values of the sample type TP. For reshaping all sample-level log2 
#' mRNASeq expression values of one cohort, see \code{\link{dn_mRNASeq_cohort}}, use \code{\link{reshape.mRNASeq_geneID}.
#' @param data data.frame of log2 mRNASeq expression values. See \code{\link{dn_mRNASeq}} for downloading these values.
#' @return nxp-matrix, n = number of patients, p = number of genes, \cr
#'         (matrix(i,j))_{i,j} = (matrix(i,j))_{i,j} =  (expression_log2)_{ij}, where \cr
#'         data$tcga_participant_barcode==barcode[i], data$gene==gene[j], \cr
#'         i=1,..,n, j=1,...,p
#' @export
#' @examples
#' gene = c("A1BG", AAAS")
#' cohort = "ESCA"
#' tcga_participant_barcode = c("TCGA-2H-A9GF", "TCGA-LN-A49M") #  TCGA patient barcodes from ESCA
#' sort_by = "gene"
#' page.Size = 250
#' obj = dn_mRNASeq(gene, cohort, tcga_participant_barcode, sort_by, page.Size)
#' obj_reshaped = reshape.mRNASeq(obj, sample_type = "TP")

reshape.mRNASeq = function(data, sample_type = "TP"){

  cat("Info: the sample types of the data:", unique(data$sample_type))
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
    matrix[,j] =  as.numeric(data[(data$gene==gene[j]), idx])
  }

  return(matrix)

}




#' Reshape sample-level log2 mRNASeq expression values of the sample type TP. Primarily this service is intended to reshape all sample-level log2 mRNASeq expression values of one cohort.
#' See \code{\link{dn_mRNASeq_cohort}}. \code{\link{gene_geneID} is a data.frame providing genes and the corresponding gene ID's.
#' @param data data.frame of log2 mRNASeq expression values. See \code{\link{dn_mRNASeq}} for downloading these values.
#' @return nxp-matrix, n = number of patients, p = number of gene ID's, \cr
#'         (matrix(i,j))_{i,j} = (matrix(i,j))_{i,j} =  (expression_log2)_{ij}, where \cr
#'         data$tcga_participant_barcode==barcode[i], data$geneID==geneID[j], \cr
#'         i=1,..,n, j=1,...,p
#' @export
#' @examples
#' cohort = "ESCA"
#' page_size = 2000
#' esca.mRNASeq = dn_mRNASeq_cohort(cohort, page.Size)
#' esca.mRNASeq_reshaped = reshape.mRNASeq_geneID(esca.mRNASeq, sample_type = "TP")
reshape.mRNASeq_geneID = function(data, sample_type = "TP"){
  
  cat("Info: the sample types of the data:", unique(data$sample_type))
  data = subset(data, data$sample_type=="TP")
  idx = which(colnames(data)=='expression_log2')
  barcode = unique(data[,"tcga_participant_barcode"])
  geneID = unique(data[,"geneID"])
  n = length(barcode)
  p = length(geneID)
  matrix = matrix(0, nrow = n, ncol = p)
  colnames(matrix) = geneID
  rownames(matrix) = barcode
  
  # fill in column by column:
  for (j in 1:p) {
    matrix[,j] =  as.numeric(data[(data$geneID==geneID[j]), idx])
  }
  
  return(matrix)
  
}

