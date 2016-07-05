
require("FirebrowseR")
require("plyr")
require("parallel")

# 1.
#'  Download all available cohorts / cancer types.
#'  @return A character vektor containing all TCGA cohort abbreviations which are relevant for the algorithms.
#'  @export
#'  @seealso \code{dn_cohorts} uses the service \code{\link{Metadata.Cohorts}}, see \code{\link{FirebrowseR}}.
dn_cohorts = function( filename=NULL ){

cohorts = Metadata.Cohorts( format = "csv", cohort = "" )
cohorts = cohorts[,1]

#  if(is.null(filename)) {
#    filename = sprintf("cohorts.Rdata");
#  }
#  save(list=c("cohorts"), file=filename);

return(cohorts)

}

# At the present: 38 cohorts
# Remark: GBMLGG is composed of GBM and LGG, COADREAD is composed of COAD and READ
#         and KIPAN is composed of KICH, KIRC and KIRP

# 2.
#' Download all available clinical data of one cohort.
#' @param cohort A character vector indicating the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @return data.frame of all patient clinical data elements of one cohort.
#' @export
#' @seealso \code{dn_clinical} uses the service \code{\link{Samples.Clinical}}, see \code{\link{FirebrowseR}}.

dn_clinical_cohort = function(cohort, filename=NULL) {

  all.Received = F
  page.Counter = 1
  page.Size = 2000
  cohort.clinical = list()

  while(all.Received == F){
    cohort.clinical[[page.Counter]] = Samples.Clinical(format = "csv", cohort = cohort, tcga_participant_barcode = "",  page = page.Counter, page_size = page.Size, sort_by = "tcga_participant_barcode")

    if(page.Counter > 1)
      colnames(cohort.clinical[[page.Counter]]) = colnames(cohort.clinical[[page.Counter-1]])

    if(nrow(cohort.clinical[[page.Counter]]) < page.Size){
      all.Received = T
    } else{
      page.Counter = page.Counter + 1
    }
  }

  cohort.clinical = do.call(rbind.fill, cohort.clinical)

 # if(is.null(filename)) {
 #  filename = sprintf("clinical_%s.rdata", toupper(cohort));
 # }

 #  save(list=c("clinical_%s", toupper(cohort)), file=filename);

  return(cohort.clinical)
}



# 3.
#' Download available clinical data of multiple / all cohorts.
#' @param cohorts A character vector indicating the cohort(s) to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @return data.frame of patient clinical data of multiple / all cohorts.
#' @export
#' @seealso \code{dn_clinical} uses the service \code{\link{dn_clinical_cohort}}. Note that for downloading all available barcodes, see \code{\link{dn_patient_barcodes}}, the object cohorts is a character vector containing all available cohort abbreviations.
dn_clinical = function(cohorts, filename=NULL) {

    clinical = list()
    page.Size = 2000

    clinical = mclapply(cohorts, dn_clinical_cohort, mc.cores = 20) # mc.scores < 32

    clinical = do.call(rbind.fill, clinical)

 #   if(is.null(filename)) {
 #     filename = sprintf("all.clinical.Rdata");
 #   }
 #  save(list=c("clinical_data"), file=filename);

   return (clinical)
}




# 4.
#' Download patient barcodes with corresponding cohort abbreviations of cohort(s) to query.
#' @param cohort A character vector indicating the cohort(s) to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @return A data.frame containing TCGA barcodes of the cohort(s) with corresponding cohort abbreviations.
#' @export
#' @seealso \code{\link{dn_clinical_all}}
dn_patient_barcodes = function(cohort, filename=NULL){

  clinical = dn_clinical(cohorts)
  barcodes = clinical[,1:2]

 #  if(is.null(filename)) {
 #    filename = sprintf("all.pats.Rdata");
 #  }
 #  save(list=c("patient_barcodes"), file=filename);

 return(barcodes)
}






# 5.
#' Download all available genes.
#' @param tcga_participant_barcode A character containing a random TCGA patient barcode. See e.g. \code{\link{dn_patient_barcodes}} for available barcodes.
#' @param page.Size Number of records per page. Usually max is 2000.
#' @return A character vector of all available gene symbols.
#' @export
#' @seealso \code{dn_gene_all} uses the service \code{\link{Analyses.CopyNumber.Genes.All}}, see \code{\link{FirebrowseR}}.


dn_gene_ID = function(tcga_participant_barcode, page.Size, filename=NULL){

  all.Received = F
  page.Counter = 1
  gene_ID = list()

  while(all.Received == F){
    gene_ID[[page.Counter]] = Analyses.CopyNumber.Genes.All(format = "csv", cohort = "", gene = "",
                                                              tcga_participant_barcode = tcga_participant_barcode, page = page.Counter,
                                                              page_size = page.Size, sort_by = "gene")
    if(page.Counter > 1)
      colnames(gene_ID[[page.Counter]]) = colnames(gene_ID[[page.Counter-1]])

    if(nrow(gene_ID[[page.Counter]]) < page.Size){
      all.Received = T
    } else{
      page.Counter = page.Counter + 1
    }
  }

  gene_ID = do.call(rbind, gene_ID)
  gene_ID = gene_ID[ ,colnames(gene_ID)=="gene"]
  idx = grep("TRUE", duplicated(gene_ID))
  gene_ID = gene_ID[-idx]

  #  if(is.null(filename)) {
  #    filename = sprintf("gene_ID.Rdata");
  #  }
  #  save(list=c("gene_IDs"), file=filename);

  return(gene_ID)
}

#' @examples
#' t = proc.time()
#' cohorts = dn_cohorts()
#' cohort = "BRCA"
#' brca.clinical = dn_clinical_cohort(cohort)
#' brca.barcodes = dn_patient_barcodes(cohort) # returns all patient barcodes of the cohort BRCA
#' clinical = dn_clinical_all("ACC", "BLCA", "BRCA") # returns clinical data of cohorts ACC, BLCA and BRCA
#' barcodes = dn_patient_barcodes(clinical) # returns all patient barcodes of cohorts ACC, BLCA and BRCA
#' tcga_participant_barcode = "TCGA-E9-A2JT" # it is a BRCA patient barcode
#' page.Size = 2000
#' gene_ID = dn_gene_ID(tcga_participant_barcode, page.Size)
#' cat(proc.time() - t)





