% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dn_mutation_SMG.R
\name{dn_mutation.SMG_cohort}
\alias{dn_mutation.SMG_cohort}
\title{Download selected columns from the MAF generated by MutSig for one cohort.}
\usage{
dn_mutation.SMG_cohort(cohort, gene_ID, page.Size, filename = NULL)
}
\arguments{
\item{cohort}{A character vector indicating the cohort to query. See \code{\link{dn_cohorts}} for available cohorts.}

\item{gene_ID}{A character vector of all available gene symbols.}

\item{page.Size}{Number of records per page. Usually max is 2000. page.Size should be chosen bigger than the number of patients of the cohort.}
}
\value{
data.frame of selected columns from the MAF generated by MutSigof for all patient barcodes of the specified cohort.
}
\description{
Download selected columns from the MAF generated by MutSig for one cohort.
}
\examples{
cohort = "ACC"
page.Size = 1000
acc.mutation.SMG = dn_mutation.SMG_cohort(cohort, page.Size)
}
\seealso{
\code{dn_mutation_cohort} uses the service \code{\link{dn_mutation.Exp}}.
}

