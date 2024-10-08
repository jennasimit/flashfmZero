% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLASHFMwithFINEMAP.R
\name{finemap}
\alias{finemap}
\title{Wrapper to run FINEMAP  (Benner et al. 2016) in R}
\usage{
finemap(GWAS, ldfile, N, maxcv = 10, fstub, FMpath)
}
\arguments{
\item{GWAS}{a data.frame with columns: "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"; this is the same as z file of FINEMAP}

\item{ldfile}{path to file that contains the SNP correlation matrix, with SNPs in the same order as in the GWAS data.frame; this is the same as the ld file of FINEMAP}

\item{N}{sample size for trait}

\item{maxcv}{maximum number of causal variants}

\item{fstub}{file stub for input/output files of FINEMAP, e.g. if fstub="DIRresults/region1", FINEMAP files of the form "DIRresults/region1.z"  will be created}

\item{FMpath}{file pathway to FINEMAP software e.g. "/software/finemap_v1.4_x86_64/finemap_v1.4_x86_64"}
}
\value{
snpPP a data.frame of top SNP models and their model PPs, as output from FINEMAP
}
\description{
Wrapper to run FINEMAP  (Benner et al. 2016) in R
}
\author{
Jenn Asimit
}
