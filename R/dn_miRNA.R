
require("FirebrowseR")
require("parallel")
require("plyr")


  

#' Download sample-level log2 miRSeq expression values. Results may be filtered by mir, cohort and barcode. At least one miR must be supplied.
#' @param mir A character vector containing miR names. At least one mir must be supplied. See \code{\link{miRNA_ID.R}} for available miR names.
#' @param cohort A character vector containing the cohort(s) to query, empty string queries all cohorts. See \code{\link{dn_cohorts}} for available cohorts.
#' @param tcga_participant_barcode A character vector containing TCGA barcodes, empty string queries all barcodes. See \code{\link{patient_barcodes}} for available barcodes. Remark that the data are NULL for barcode(s) which isn´t (aren´t) barcode(s) of the specified cohort.
#' @param sort_by A character indicating the column which is used for sorting. The data can be sorted by tcga_participant_barcode, cohort, tool, mir and sample_type.
#' @param page.Size Number of records per page. Usually max is 2000.
#' @return data.frame of log2 miRSeq expression values.
#' @export
#' @seealso \code{dn_miRSeq} uses the service \code{\link{Samples.miRSeq}}, see \code{\link{FirebrowseR}}.
#' @examples
#' mir = miRNA_ID[1:10]
#' cohort = "BLCA"
#' tcga_participant_barcode = "TCGA-ZF-AA53"  # it is a BLCA patient barcode
#' page.Size = 250
#' sort_by =  "tcga_participant_barcode"
#' obj = dn_miRSeq(mir, cohort, tcga_participant_barcode, page.Size, sort_by)
dn_miRSeq = function(mir, cohort, tcga_participant_barcode, page.Size, sort_by, filename=NULL){

  all.Found = F
  page.Counter = 1
  mir.Exp = list()

  while(all.Found == F){
    tmp = Samples.miRSeq(format = "csv", mir = mir, cohort = cohort,
                         tcga_participant_barcode = tcga_participant_barcode, tool = "miRseq_Mature_Preprocess", sample_type = "",
                         page = page.Counter, page_size = page.Size, sort_by = sort_by)

    if( is.null(tmp)==TRUE) { tmp = NULL; break; }

    mir.Exp[[page.Counter]] = tmp

    if(page.Counter > 1){
      colnames(mir.Exp[[page.Counter]]) = colnames(mir.Exp[[page.Counter-1]])
    }

    if(nrow(mir.Exp[[page.Counter]]) < page.Size){
      all.Found = T
    } else{
      page.Counter = page.Counter + 1
    }
  }

  if(length(mir.Exp)<1) {
    mir.Exp = NULL
  } else{
    mir.Exp = do.call(rbind.fill, mir.Exp)
  }

  return(mir.Exp)

}




#' Download all available sample-level log2 miRSeq expression values of one cohort.
#' @param cohort A character vector indicating the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @param page.Size Number of records per page. Usually max is 2000.
#' @return data.frame of sample-level log2 miRSeq expression values of the specified cohort.
#' @export
#' @seealso \code{dn_miRSeq_cohort} uses the service \code{\link{dn_miRSeq}}.
#' @examples
#' cohort = "BLCA"
#' page.Size = 2000
#' blca.miRSeq = dn_miRSeq_cohort(cohort, page.Size)
dn_miRSeq_cohort = function(cohort, page.Size, filename=NULL){

  cohort.miRSeq = list()

  cohort.miRSeq = mclapply(miRNA_ID, dn_miRSeq, cohort, "", page.Size, "mir", mc.cores=detectCores()/2)

  if (length(cohort.miRSeq) < 1){
    cohort.miRSeq = NULL
  } else{
    cohort.miRSeq= do.call(rbind.fill, cohort.miRSeq)
  }

  return(cohort.miRSeq)

}



