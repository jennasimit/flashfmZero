% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/No-Cor-flashfm-v20240222_update.R
\name{calcQ12}
\alias{calcQ12}
\title{internal function for calcAdjPP for a pair of traits}
\usage{
calcQ12(i, j, T1, T2, SS, PP, tau, nsnpspermodel, kappa)
}
\arguments{
\item{i}{model index for trait 1}

\item{j}{model index for trait 2}

\item{T1}{index of trait 1}

\item{T2}{index of trait 2}

\item{SS}{list consisting of lists of model configuration SNPs for each trait}

\item{PP}{list consisting of lists of model PP for each trait}

\item{tau}{matrix of adjustment terms}

\item{nsnpspermodel}{list of number of SNPs per model for each model in STR}

\item{kappa}{single value of sharing parameter kappa}
}
\description{
internal function for calcAdjPP for a pair of traits
}
\author{
Jenn Asimit
}
