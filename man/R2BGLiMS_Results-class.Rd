% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classes.R
\docType{class}
\name{R2BGLiMS_Results-class}
\alias{R2BGLiMS_Results-class}
\alias{R2BGLiMS_Results,R2BGLiMS_Results-class}
\title{The R2BGLiMS_Results class}
\description{
Container for results from running R2BGLiMS as generated
by the function R2BGLiMS.
}
\section{Slots}{

\describe{
\item{\code{likelihood}}{Likelihood type of model. See R2BGLiMS for the options.}

\item{\code{posterior.summary.table}}{Posterior summaries of all parameters}

\item{\code{enumerate.up.to.dim}}{If posterior inference was made by exhaustively enumerating and assessing models one by one, modles were considered
up to this dimension.}

\item{\code{n.iterations}}{Number of iterations which the RJMCMC was run for.}

\item{\code{thin}}{Ith iterations which were saved.}

\item{\code{model.space.prior}}{List defining the model space prior. See R2BGLiMS.}

\item{\code{beta.prior.partitions}}{List defining the covariate effect prior partitions. See R2BGLiMS.}

\item{\code{confounders}}{Vectors of variables fixed in the model, and excluded from model selection.}

\item{\code{run.times}}{A list containing run times broken down into different processes.}

\item{\code{n.covariate.blocks.for.jam}}{The number of partitioned LD blocks used for JAM.}

\item{\code{bglims.arguments}}{The arguments passed to the Java BGLiMS function. This is a named list - the different arguments are a range of datatypes.}

\item{\code{mcmc.output}}{The Reversible Jump MCMC output from BGLiMS. Columns are parameters, rows are iterations.}
}}

\author{
Paul J. Newcombe \email{paul.newcombe@mrc-bsu.cam.ac.uk}
}
