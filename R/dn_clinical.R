
require("FirebrowseR")
require("plyr") 
require("parallel")



#' Download all available cohorts / cancer types.
#' @return A character vektor containing all TCGA cohort abbreviations which are relevant for the algorithms.
#' @export
#' @seealso \code{dn_cohorts} uses the service \code{\link{Metadata.Cohorts}}, see \code{\link{FirebrowseR}}.
#' @examples
#' cohorts = dn_cohorts()
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
#' @seealso \code{dn_clinical_one} uses the service \code{\link{Samples.Clinical}}, see \code{\link{FirebrowseR}}.
#' @examples
#' cohort = "BRCA"
#' brca.clinical = dn_clinical_one(cohort)
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
#' @seealso \code{dn_clinical} uses the service \code{\link{dn_clinical_one}}. Note that for downloading all available barcodes, see \code{\link{patient_barcodes}}, the object cohorts is a character vector containing all available cohort abbreviations.
#' @examples
#' cohorts = c("ACC", "BLCA", "COAD")
#' clinical = dn_clinical(cohorts)
dn_clinical = function(cohorts, filename=NULL) {

    clinical = list()
    page.Size = 2000

    clinical = mclapply(cohorts, dn_clinical_one, mc.cores = detectCores()/2)

    clinical = do.call(rbind.fill, clinical)

 #   if(is.null(filename)) {
 #     filename = sprintf("clinical.Rdata");
 #   }
 #  save(list=c("clinical_data"), file=filename);

   return (clinical)  
}
 
  


#' Extract TCGA patient barcodes from clinical data.
#' @param clinical A data frame containing clinical data. See \code{\link{dn_clinical_one}} and \code{\link{dn_clinical}}.
#' @return A character vector containing TCGA barcodes.
#' @export
#' @examples
#' cohort = "BRCA"
#' brca.clinical = dn_clinical_one(cohort)
#' brca.barcodes = patient_barcodes(clinical) # returns all patient barcodes of the cohort BRCA
patient_barcodes = function(clinical){

  barcodes = clinical[,"tcga_participant_barcode"]

  return(barcodes)
}


