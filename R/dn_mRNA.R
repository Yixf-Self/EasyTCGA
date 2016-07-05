require("FirebrowseR")
require("parallel")
# source("dn_clinical.R")


#--------------------------------------------------------------------------------
# This algorithm doens´t work because of the updated version
# 1.
#' Download all mRNA IDs.
#' @param tcga_participant_barcode A character containing a random TCGA patient barcode, e.g. the barcode "TCGA-E9-A2JT" can be used.
#' It may happen that the data are NULL, the reasons can be:
#'     - There are no data for the specified barcode
#'     - For some cancer types (cohorts) there are no data at all -> try a tcga_participant_barcode of another cancer type
#'       Problem : without mRNA IDs there is no possibility to find out which cancers don´t have mRNARSeq expression values
#'       The good news: the output is the same for any adequate barcode
#'     - Lack of connection to the server -> try the same barcode again or choose another barcode
#'     -> please check if mRNA_ID, the output of \code{dn_mRNA_ID}, is NULL for the specified barcode before using the service
#'        \code{\link{dn_mRNASeq_cohort}}, otherwise the output of \code{\link{dn_mRNASeq_cohort}} is NULL as well
#' @param page.Size Number of records per page. Usually max is 2000.
#' @return A character vector containing all available mRNA IDs.
#' @export
#' @seealso \code{dn_mRNA_ID} uses the service \code{\link{Samples.mRNASeq}}, see \code{\link{FirebrowseR}}.

dn_mRNA_ID = function(tcga_participant_barcode, page.Size, filename=NULL){

  all.Found = F
  page.Counter = 1
  mRNA_ID = list()

  while(all.Found == F){
    tmp = Samples.mRNASeq(format = "csv", gene = "", cohort = "", tcga_participant_barcode = tcga_participant_barcode,
                          sample_type = "", protocol = "RSEM", page = page.Counter, page_size = page.Size, sort_by = "gene")

    if( is.null(tmp)==TRUE) { tmp = NULL; break; }

    mRNA_ID[[page.Counter]] = tmp


    if(page.Counter > 1){
      colnames( mRNA_ID[[page.Counter]]) = colnames( mRNA_ID[[page.Counter-1]])
    }

    if(nrow( mRNA_ID[[page.Counter]]) < page.Size){
      all.Found = T
    } else{
      page.Counter = page.Counter + 1
    }

  }

  mRNA_ID = do.call(rbind.fill, mRNA_ID)

  # The second column of the downloaded data contains all mRNAs:
  mRNA_ID = mRNA_ID[,2]

  # Eliminate duplicated mRNA IDs:
  idx1 = grep("TRUE", duplicated(mRNA_ID))
  idx2 = grep("TRUE", mRNA_ID=="None")
  idx = union(idx1, idx2)
  mRNA_ID = mRNA_ID[-idx]


  #  if(is.null(filename)) {
  #    filename = sprintf("miRNA_IDs.rdata");
  #  }

  # output is stored in miRNA_ID <- miRNA_ID
  #  save(list=c("miRNA_IDs"), file=filename);

  return(mRNA_ID)
}

# save(file="mRNA_ID.RData", list="mRNA_ID")
#----------------------------------------------------------------------------------


# 2.
#' Download sample-level log2 mRNASeq expression values. Results may be filtered by gene, tcga_participant_barcode and cohort. At least one gene must be supplied.
#' @param tcga_participant_barcode A character vector containing TCGA barcodes, empty string queries all barcodes. See e.g. \code{\link{dn_patient_barcodes}} for available barcodes. Remark that the data are NULL for barcode(s) which isn´t (aren´t) barcode(s) of the specified cohort.
#' @param cohort A character vector indicating cohort(s) to query, empty string queries all cohorts. See \code{\link{dn_cohorts}} for available cohorts.
#' @param gene A character vector of gene symbols. At least one gene must be supplied. See \code{\link{mRNA_ID.R}} for available genes.
#' @param sort_by A character indicating the column which is used for sorting. The data can be sorted by tcga_participant_barcode, cohort, gene, protocol and sample_type.
#' @param page.Size Number of records per page. Usually max is 2000. For cancers with small number of patients use a small page.Size, e.g. page.Size = 250. # If there is an error "Error in mRNA.Exp[[page.Counter]]: subscript out of bounds", run the algorithm again or try a smaller page.Size.
#' @return data.frame of log2 mRNASeq expression values.
#' @export
#' @seealso  \code{dn_mRNASeq} uses the service \code{\link{Samples.mRNASeq}}, see \code{\link{FirebrowseR}}.

dn_mRNASeq = function(gene, cohort, tcga_participant_barcode, sort_by, page.Size, filename=NULL) {

  all.Found = F
  page.Counter = 1
  mRNA.Exp = list()
  # cat("gene = ", gene, "cohort= ", cohort, " tcga_participant_barcode= ", tcga_participant_barcode, " sort_by= ", sort_by, "page.Size= ", page.Size, \n")

#  cat(".")
#  flush.console()

  while(all.Found == F){
    tmp = Samples.mRNASeq(format = "csv", gene = gene, cohort = cohort,
                                               tcga_participant_barcode = tcga_participant_barcode, sample_type = "", protocol = "RSEM",
                                               page = page.Counter, page_size = page.Size, sort_by = sort_by)

    if( is.null(tmp)==TRUE) {tmp = NULL; break;}

    mRNA.Exp[[page.Counter]] = tmp

    if(page.Counter > 1){
      colnames(mRNA.Exp[[page.Counter]]) = colnames(mRNA.Exp[[page.Counter-1]])
    }

    if(nrow(mRNA.Exp[[page.Counter]]) < page.Size)
      all.Found = T
    else
      page.Counter = page.Counter + 1
  }

  if(length(mRNA.Exp) < 1) {
    mRNA.Exp = NULL
  } else {
    mRNA.Exp = do.call(rbind.fill, mRNA.Exp)
  }
 # cat("\n")

#  if(is.null(filename)) {
#    filename = sprintf("mRNA.Exp.rdata");
#  }

 # cat("o")
#  flush.console()

  # output is stored in clinical <- cancerType.Pats
  # save(list=c("mRNA.Exp", "cohort", "gene", "tcga_participant_barcode"), file=filename);

  return(mRNA.Exp)
}



# 3.
#' Download all available sample-level log2 mRNASeq expression values of one cohort.
#' @param mRNA_ID A character vector containing gene symbols. See \code{\link{mRNA_ID.R}} for available gene symbols.
#' @param cohort A character vector containing the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @param page.Size Number of records per page. Usually max is 2000. For cancers with small number of patients use a small page.Size, e.g. page.Size = 250.
#' @return data.frame of sample-level log2 mRNASeq expression values of the specified cohort.
#' @export
#' @seealso \code{dn_mRNASeq_cohort} uses the service \code{\link{dn_mRNASeq}}.


dn_mRNASeq_cohort = function(cohort, page.Size, filename=NULL){

  # cohort.clinical = dn_clinical(cohort)
  # cohort.mRNASeq = mclapply(mRNA_ID, dn_mRNASeq, cohort, "", "tcga_participant_barcode", page.Size, mc.cores=32)    # <32
  cohort.mRNASeq = mclapply(mRNA_ID, dn_mRNASeq, cohort, "", "gene", page.Size, mc.cores=32)
  # print(mRNA.cohort)

  if(length(cohort.mRNASeq)<1){
    cohort.mRNASeq = NULL
  } else{
    cohort.mRNASeq = do.call(rbind, cohort.mRNASeq)
  }

  return(cohort.mRNASeq)

}







#' @examples
#' t = proc.time()
#'
#' gene = "AAAS" # "=mRNA_ID[10]"
#' cohort = "ESCA"
#' tcga_participant_barcode = c("TCGA-2H-A9GF", "TCGA-LN-A49M") # esca patient barcodes
#' sort_by = "gene"
#' page.Size1 = 250
#' obj = dn_mRNASeq(gene, cohort, tcga_participant_barcode, sort_by, page.Size2)
#' page_size2 = 2000
#' esca.mRNASeq = dn_mRNASeq_cohort(cohort, page.Size2)
#'
#' cat(proc.time() - t)

