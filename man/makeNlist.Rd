% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flashfm.R
\name{makeNlist}
\alias{makeNlist}
\title{Sample size information needed for flashfm}
\usage{
makeNlist(Nall, y = NULL, Nsame = NULL)
}
\arguments{
\item{Nall}{vector of sample sizes for each trait}

\item{y}{if avaialable, matrix of trait measurements or indicators of non-NA trait measurements (columns are traits);
used to get joint sample counts; default is NULL and if not provided an approximation is used based on vector of trait sample sizes}

\item{Nsame}{a single sample size that is the same for all traits}
}
\value{
list of 4 components: N = number of individuals with all traits measured;
Nqq=matrix of all pair-wise counts of number of individuals with both traits in a pair measured;
Nq3 = vector of counts of number of individuals with three traits measured; all triples considered; NULL if M < 4
Nq4 = vector of counts of number of individuals with four traits measured; all quadruples considered; NULL if M < 5
}
\description{
Sample size information needed for flashfm
}
\author{
Jenn Asimit
}
