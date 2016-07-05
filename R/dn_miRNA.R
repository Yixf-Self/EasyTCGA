require("FirebrowseR")
require("parallel")
# source("dn_clinical.R")

#--------------------------------------------------------------
# This algorithm doens´t work because of the updated version
#' Download all miR names.
#' @param tcga_participant_barcode A character containing a random TCGA patient barcode, the barcode "TCGA-BH-A1EO" can be used. See \code{\link{dn_patient_barcodes}} for available barcodes.
#' @return A character vector containing all available miR names.
#' @export
#' @seealso \code{dn_miRNA_ID} uses the service \code{\link{Samples.miRSeq}}, see \code{\link{FirebrowseR}}.
#'   It may happen that the data are NULL, the reasons can be:
#'     - There are no data for the specified barcode
#'     - For some cancer types (cohorts) there are no data at all
#'       -> try a tcga_participant_barcode of another cancer type
#'       Problem : without miRNA IDs there is no possibility to find out which cancers don´t have miRSeq expression values
#'       The good news: the output is the same for any adequate barcode
#'     - Lack of connection to the server -> try the same barcode again or choose another barcode


dn_miRNA_ID = function(tcga_participant_barcode, filename=NULL){

  all.Found = F
  page.Counter = 1
  miRNA_ID = list()
  page.Size = 2000

  while(all.Found == F){
    tmp = Samples.miRSeq(format = "csv", mir = "", cohort = "",
                         tcga_participant_barcode = tcga_participant_barcode, tool = "miRseq_Mature_Preprocess",
                         sample_type = "", page = page.Counter, page_size = page.Size, sort_by = "mir")

    # This if-statement and the following command "miRNA_ID[[page.Counter]] = tmp" as an intermediate step is necessary for
    # the page.Counter = 1: if the data are NULL for the specified barcode then using "miRNA_ID[[page.Counter]] = Samples.miRSeq(...)"
    # directly causes the error "Error in miRNA_ID[[page.Counter]]: subscript out of bounds" while running the if-else-statement below
    #   if(nrow( miRNA_ID[[page.Counter]]) < page.Size){
    #    all.Found = T
    #   } else{
    #    page.Counter = page.Counter + 1
    #   }
    if( is.null(tmp)==TRUE) {
      tmp = NULL; break;
    }

    miRNA_ID[[page.Counter]] = tmp

    if(page.Counter > 1){
      colnames( miRNA_ID[[page.Counter]]) = colnames( miRNA_ID[[page.Counter-1]])
    }

    if(nrow( miRNA_ID[[page.Counter]]) < page.Size){
      all.Found = T
    } else{
      page.Counter = page.Counter + 1
    }

  }

  if(length(miRNA_ID)<1){
    miRNA_ID = NULL
  } else{
    miRNA_ID = do.call(rbind.fill, miRNA_ID)
  }

  # The second column of the downloaded data contains all miRNAs
  miRNA_ID = miRNA_ID[,2]
  # In some cases the same IDs are spelled differently
  miRNA_ID = tolower(miRNA_ID)
  idx = grep("TRUE", duplicated(miRNA_ID))
  miRNA_ID = miRNA_ID[-idx]

  #  if(is.null(filename)) {
  #    filename = sprintf("miRNA_IDs.rdata");
  #  }
  #  save(list=c("miRNA_IDs"), file=filename);

  return(miRNA_ID)
}

# Save the IDs in advance:
# save(file="miRNA_ID.RData", list="miRNA_ID")
#-------------------------------------------------------------------------------------------



# 2.
#' Download sample-level log2 miRSeq expression values. Results may be filtered by mir, cohort and barcode. At least one miR must be supplied.
#' @param tcga_participant_barcode A character vector containing TCGA barcodes, empty string queries all barcodes. See \code{\link{dn_patient_barcodes}} for available barcodes. Remark that the data are NULL for barcode(s) which isn´t (aren´t) barcode(s) of the specified cohort.
#' @param cohort A character vector containing the cohort(s) to query, empty string queries all cohorts. See \code{\link{dn_cohorts}} for available cohorts.
#' @param mir A character vector containing miR names, empty string queries all miR names. See \code{\link{miRNA_ID.R}} for available miR names.
#' @param sort_by A character indicating the column which is used for sorting. The data can be sorted by tcga_participant_barcode, cohort, tool, mir and sample_type.
#' @param page.Size Number of records per page. Usually max is 2000.
#' @return data.frame of log2 miRSeq expression values.
#' @export
#' @seealso  \code{dn_miRSeq} uses the service \code{\link{Samples.miRSeq}}, see \code{\link{FirebrowseR}}.

#dn_miRSeq = function(tcga_participant_barcode, cohort, mir, sort_by, page.Size, filename=NULL){
dn_miRSeq = function(mir, cohort, tcga_participant_barcode, sort_by, page.Size, filename=NULL){
  all.Found = F
  page.Counter = 1
  mir.Exp = list()

  while(all.Found == F){
    tmp = Samples.miRSeq(format = "csv", mir = mir, cohort = cohort,
                         tcga_participant_barcode = tcga_participant_barcode, tool = "miRseq_Mature_Preprocess", sample_type = "",
                         page = page.Counter, page_size = page.Size, sort_by = sort_by)

    # This if-statement and the following command "miRNA.Exp[[page.Counter]] = tmp" as an intermediate step is necessary for
    # the page.Counter = 1: if the data are NULL for the specified barcode and/or gene then using "miRNA.Exp[[page.Counter]] = Samples.miRSeq(...)"
    # directly causes the error "Error in miRNA.Exp[[page.Counter]]: subscript out of bounds" while running the if-else-statement below
    #   if(nrow( miRNA.Exp[[page.Counter]]) < page.Size){
    #    all.Found = T
    #   } else{
    #    page.Counter = page.Counter + 1
    #   }

    #
    #cat(tmp)
    #

    if( is.null(tmp)==TRUE) { tmp = NULL; break; }

    mir.Exp[[page.Counter]] = tmp

    if(page.Counter > 1){
      colnames(mir.Exp[[page.Counter]]) = colnames(mir.Exp[[page.Counter-1]])
    }

    #
    #nrow(mir.Exp[[page.Counter]])
    #

    if(nrow(mir.Exp[[page.Counter]]) < page.Size){
      all.Found = T
    } else{
      page.Counter = page.Counter + 1
    }
  }

  #
  #length(mir.Exp)
  #

  if(length(mir.Exp)<1) {
    mir.Exp = NULL
  } else{
    mir.Exp = do.call(rbind.fill, mir.Exp)
  }

  #if(is.null(filename)) {
  #  filename = sprintf("miRSeq.rdata");
  #}

  # save(list=c("miRSeq", "cohort", "tcga_participant_barcode", "miRNA_ID"), file=filename);
  # cat(".")

  return(mir.Exp)

}




# 3.
#' Download all available sample-level log2 miRSeq expression values of one cohort.
#' @param cohort A character vector indicating the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @param page.Size Number of records per page. Usually max is 2000. # If there is an error "Error in miRNA.cohort[[page.Counter]]: subscript out of bounds", try a smaller page.Size.
#' @return data.frame of sample-level log2 miRSeq expression values of the specified cohort.
#' @export
#' @seealso \code{dn_miRNA_cohort} uses the service \code{\link{dn_miRSeq}}.
#'
dn_miRSeq_cohort = function(cohort, page.Size, filename=NULL){

  cohort.miRSeq = list()
  # cohort.clinical = dn_clinical(cohort)
  # dn_miRSeq = function(tcga_participant_barcode, cohort, mir, sort_by, page.Size, filename=NULL)
  # miRNA.cohort = mclapply(cohort.clinical[,1], dn_miRSeq, cohort, "", "tcga_participant_barcode", page.Size, mc.cores=20)
  cohort.miRSeq = mclapply(miRNA_ID, dn_miRSeq, cohort, "", "mir", page.Size, mc.cores=20)

  if (length(cohort.miRSeq)<1){
    cohort.miRSeq = NULL
  } else{
    cohort.miRSeq= do.call(rbind.fill, cohort.miRSeq)
  }

  #  if(is.null(filename)) {
  #    filename = sprintf("%s.miRSeq.rdata", toupper(cohort));
  #  }

  #  save(list=c("", "cohort"), file=filename);

  return(cohort.miRSeq)

}





# str = sprintf("Mut_%s.rdata", cohort)     # %d integer %f float ...
# mutation = Mut.cohort
# save(file=str, list=c("mutation", "cohort"))




#'@examples
#' t = proc.time()
#'
#' cohort = "BLCA"
#' tcga_participant_barcode = "TCGA-ZF-AA53"  # it is a BLCA patient barcode
#' page.Size1 = 250
#' sort_by =  "tcga_participant_barcode"
#' mir = miRNA_ID[1:10]
#' obj = dn_miRSeq(mir, cohort, tcga_participant_barcode, sort_by, page.Size1)
#' page.Size2 = 2000
#' blca.miRSeq = dn_miRSeq_cohort(cohort, page.Size2)
#'
#' cat(proc.time() - t)





