% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrettyResultsTable.R
\name{PrettyResultsTable}
\alias{PrettyResultsTable}
\title{A pretty summary results table for reports. NOTE: This function outputs a table formatted with character strings.
A numeric representation of the results are stored in the slot 'posterior.summary.table'.}
\usage{
PrettyResultsTable(
  results,
  round.digits.betas = 2,
  round.digits.postprob = 2,
  round.digits.bf = 1,
  normalised.sds = NULL
)
}
\arguments{
\item{results}{A R2BGLiMS R2BGLiMS_Results-class object frum running the R2BGLiMS R2BGLiMS command.}

\item{round.digits.betas}{Number of decimal places to include for effect estimates. (Default is 2)}

\item{round.digits.postprob}{Number of decimal places to include for posterior probabilities. (Default is 2)}

\item{round.digits.bf}{Number of decimal places to include for Bayes factors. (Default is 1)}

\item{normalised.sds}{If covariates were normalised by their standard deviation, provide a named
vector of the standard deviations on the original scale using this argument. Effects in the resulting
table will then be interpretable according to unit changes on each covariate's original scale, rather
than according to standard deviation changes.}
}
\value{
A nice summary results table.
}
\description{
Creates a `pretty' summary results table from a Reversible Jump results object
}
\author{
Paul Newcombe
}
