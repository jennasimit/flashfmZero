% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flashfm.R
\name{allC12}
\alias{allC12}
\title{internal function for calcAdjPP for that gives list of covariance matrix of residuals for all trait pairs}
\usage{
allC12(M, nummods, beta, SSy, Sxy, xcovo, Mx, Nqq)
}
\arguments{
\item{M}{number of traits}

\item{nummods}{list where component i is the number of models for trait i}

\item{beta}{list of joint beta estimates for each trait}

\item{SSy}{matrix of trait cross-products}

\item{Sxy}{matrix with each column being the cross-product between SNPs and a trait}

\item{xcovo}{covariance matrix of c("one"=1,X)}

\item{Mx}{vector of SNP means}

\item{Nqq}{has  \code{N[i,j]} = no. with  both trait i and trait j measured}
}
\description{
internal function for calcAdjPP for that gives list of covariance matrix of residuals for all trait pairs
}
\author{
Jenn Asimit
}
