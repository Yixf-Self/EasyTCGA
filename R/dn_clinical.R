
require("FirebrowseR")
require("plyr")
require("parallel")

#' Download all available cohorts / cancer types.
#' @return A character vektor containing all TCGA cohort abbreviations which are relevant for the algorithms.
#' @export
#' @seealso \code{dn_cohorts} uses the service \code{\link{Metadata.Cohorts}}, see \code{\link{FirebrowseR}}.
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

#' Download all available clinical data of one cohort.
#' @param cohort A character vector indicating the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @return data.frame of all patient clinical data elements of one cohort.
#' @export
#' @seealso \code{dn_clinical} uses the service \code{\link{Samples.Clinical}}, see \code{\link{FirebrowseR}}.
dn_clinical_one = function(cohort, filename=NULL) {

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



#' Download available clinical data of multiple / all cohorts.
#' @param cohorts A character vector indicating the cohort(s) to query. See \code{\link{dn_cohorts}} for available cohorts.
#' @return data.frame of patient clinical data of multiple / all cohorts.
#' @export
#' @seealso \code{dn_clinical} uses the service \code{\link{dn_clinical_cohort}}. Note that for downloading all available barcodes, see \code{\link{dn_patient_barcodes}}, the object cohorts is a character vector containing all available cohort abbreviations.
dn_clinical = function(cohorts, filename=NULL) {

    clinical = list()
    page.Size = 2000

    clinical = mclapply(cohorts, dn_clinical_one, mc.cores = 20) # mc.scores < 32

    clinical = do.call(rbind.fill, clinical)

 #   if(is.null(filename)) {
 #     filename = sprintf("all.clinical.Rdata");
 #   }
 #  save(list=c("clinical_data"), file=filename);

   return (clinical)
}



#' Extract TCGA patient barcodes from clinical data.
#' @param clinical A data frame containing clinical data. See \code{\link{dn_clinical}}.
#' @return A character vector containing TCGA barcodes.
#' @export
patient_barcodes = function(clinical){

  barcodes = clinical[,"tcga_participant_barcode"]

 return(barcodes)
}



#' Download all available genes.
#' @param tcga_participant_barcode A character containing a random TCGA patient barcode. See e.g. \code{\link{dn_patient_barcodes}} for available barcodes.
#' @param page.Size Number of records per page. Usually max is 2000.
#' @return A character vector of all available gene symbols.
#' @export
#' @seealso \code{dn_gene_all} uses the service \code{\link{Analyses.CopyNumber.Genes.All}}, see \code{\link{FirebrowseR}}.
#' @examples
#' cohort = "BRCA"
#' clinical = dn_clinical(cohort)
#' barcodes = patient_barcodes(clinical) # returns all patient barcodes of the cohort BRCA
#' gene_IDs = dn_gene_ID(barcodes[1], page.Size=2000)
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





