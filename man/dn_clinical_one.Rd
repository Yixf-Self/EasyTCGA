% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dn_clinical.R
\name{dn_clinical_one}
\alias{dn_clinical_one}
\title{Download all available clinical data of one cohort.}
\usage{
dn_clinical_one(cohort, filename = NULL)
}
\arguments{
\item{cohort}{A character vector indicating the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.}
}
\value{
data.frame of all patient clinical data elements of one cohort.
}
\description{
Download all available clinical data of one cohort.
}
\examples{
cohort = "BRCA"
brca.clinical = dn_clinical_one(cohort)
}
\seealso{
\code{dn_clinical_one} uses the service \code{\link{Samples.Clinical}}, see \code{\link{FirebrowseR}}.
}

