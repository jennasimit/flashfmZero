% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JAMupdate.R
\name{JAMexpandedCor.multi2}
\alias{JAMexpandedCor.multi2}
\title{Expanded version of JAM (a single-trait fine-mapping approach) that first runs on thinned SNPs and then expands models on tag SNPs; this can run independently on multiple traits
This version is more stable than JAMexpandedCor.multi, but slower, so is run only if JAMexpandedCor.multi fails}
\usage{
JAMexpandedCor.multi2(
  beta1,
  corX,
  raf,
  ybar,
  Vy,
  N,
  maxcv = 10,
  r2 = 0.99,
  save.path
)
}
\arguments{
\item{beta1}{list where each component is a named vector of of single SNP effect estimates for a trait; one vector for each trait}

\item{corX}{genotype correlation matrix (reference or from sample)}

\item{raf}{named vector of reference allele frequencies; the name of each allele frequency is the SNP ID and MUST be in same SNP order as in corX}

\item{ybar}{vector of trait means; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set ybar as 0-vector}

\item{Vy}{vector of trait variances; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set Vy as 1-vector}

\item{N}{vector of sample sizes for each trait; recommended to give effective sample sizes using GWAS summary statistics in Neff function}

\item{maxcv}{maximum number of causal variants}

\item{r2}{r.squared threshold for thinning SNPs before JAM and finding tag SNPs}

\item{save.path}{path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").}
}
\value{
list with 3 components: SM a list of snpmod objects giving fine-mapping results for each trait; mbeta a list of joint effect estimates for each trait; nsnps number of SNPs
}
\description{
Expanded version of JAM (a single-trait fine-mapping approach) that first runs on thinned SNPs and then expands models on tag SNPs; this can run independently on multiple traits
This version is more stable than JAMexpandedCor.multi, but slower, so is run only if JAMexpandedCor.multi fails
}
\author{
Jenn Asimit
}
