
require("FirebrowseR")
require("parallel")
require("plyr")




#' Download sample-level log2 mRNASeq expression values. Results may be filtered by gene, tcga_participant_barcode and cohort.
#' @param gene A character vector of gene symbols, empty string or gene = mRNA_ID queries all genes. At least one gene must be supplied. See \code{\link{mRNA_ID.R}} for available genes.
#' @param cohort A character vector indicating cohort(s) to query, empty string queries all cohorts. See \code{\link{dn_cohorts}} for available cohorts.
#' @param tcga_participant_barcode A character vector containing TCGA barcodes, empty string queries all barcodes. See e.g. \code{\link{patient_barcodes}} for available barcodes. Remark that the data can be NULL for barcode(s) which isn´t (aren´t) barcode(s) of the specified cohort.
#' @param page.Size Number of records per page. Usually max is 2000. 
#' @param sort_by A character indicating the column which is used for sorting. The data can be sorted by tcga_participant_barcode, cohort, gene, protocol and sample_type.
#' @return data.frame of log2 mRNASeq expression values.
#' @export
#' @seealso  \code{dn_mRNASeq} uses the service \code{\link{Samples.mRNASeq}}, see \code{\link{FirebrowseR}}.
#' @examples
#' gene = "AAAS" 
#' cohort = "ESCA"
#' tcga_participant_barcode = c("TCGA-2H-A9GF", "TCGA-LN-A49M") #  TCGA patient barcodes from ESCA
#' sort_by = "gene"
#' page.Size = 250
#' obj = dn_mRNASeq(gene, cohort, tcga_participant_barcode, sort_by, page.Size2)
dn_mRNASeq = function(gene, cohort, tcga_participant_barcode, page.Size, sort_by, filename=NULL) {

  all.Found = F
  page.Counter = 1
  mRNA.Exp = list()

  while(all.Found == F){
  #  tmp = Samples.mRNASeq(format = "csv", gene = gene, cohort = cohort,
  #                       tcga_participant_barcode = tcga_participant_barcode, sample_type = "", protocol = "RSEM",
  #                        page = page.Counter, page_size = page.Size, sort_by = sort_by)

    tmp = tryCatch( Samples.mRNASeq(format = "csv", gene = gene, cohort = cohort,
                                    tcga_participant_barcode = tcga_participant_barcode, 
                                    sample_type = "", protocol = "RSEM", page = page.Counter, 
                                    page_size = page.Size, sort_by = sort_by) , error=function(e) NULL)
    
    if(is.null(tmp)==TRUE || length(tmp)<1) { tmp = NULL; mRNA.Exp = as.data.frame(tmp); break; }

    mRNA.Exp[[page.Counter]] = tmp

    if(page.Counter > 1){
      colnames(mRNA.Exp[[page.Counter]]) = colnames(mRNA.Exp[[page.Counter-1]])
    }

    if(nrow(mRNA.Exp[[page.Counter]]) < page.Size){
      all.Found = T
    } else{
      page.Counter = page.Counter + 1
    }
  }

  if(is.null(mRNA.Exp)==TRUE || length(mRNA.Exp) < 1) {
    mRNA.Exp = NULL
    mRNA.Exp = as.data.frame(mRNA.Exp)
  } else{
    mRNA.Exp = do.call(rbind.fill, mRNA.Exp)
  }

  return(mRNA.Exp)

}




#' Download all available sample-level log2 mRNASeq expression values of one cohort.
#' @param mRNA_ID A character vector containing gene symbols. See \code{\link{mRNA_ID.R}} for available gene symbols.
#' @param cohort A character vector containing the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @param page.Size Number of records per page. Usually max is 2000. page.Size should be chosen bigger than the number of patients of the cohort.
#' @return data.frame of sample-level log2 mRNASeq expression values of the specified cohort.
#' @export
#' @seealso \code{dn_mRNASeq_cohort} uses the service \code{\link{dn_mRNASeq}}.
#' @examples
#' cohort = "ESCA"
#' page_size = 2000
#' esca.mRNASeq = dn_mRNASeq_cohort(cohort, page.Size2)
dn_mRNASeq_cohort = function(cohort, page.Size, filename=NULL){
  
  cohort.mRNASeq = list()

  cohort.mRNASeq = mclapply(mRNA_ID, dn_mRNASeq, cohort, "", page.Size, "gene", mc.cores=detectCores()/5)

  if(is.null(cohort.mRNASeq)==TRUE || length(cohort.mRNASeq)<1){
    cohort.mRNASeq = NULL
  } else{
    cohort.mRNASeq = do.call(rbind, cohort.mRNASeq)
  }

  return(cohort.mRNASeq)

}


