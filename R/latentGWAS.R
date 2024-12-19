#' @title GWAS for latent factors at one variant using observed trait GWAS summary statistics
#' @param obsGWAS list of observed GWAS data.frames, containing the effect estimates and standard errors and snp names, at a minimum
#' @param covY trait covariance matrix (same as used in factor analysis)
#' @param L PxK matrix of factor loadings generated from a factor analysis of P observed traits and K latent factors (K<P)
#' @param beta_colname text column name for the effect estimates in each obsGWAS data.frame; default "BETA"
#' @param se_colname text column name for the effect estimate standard errors in each obsGWAS data.frame; default "SE"
#' @param snpID_colname text column name for the effect estimates in each obsGWAS data.frame; default "SNP"
#' @param EAF vector of effect allele frequencies in same order as variants in obsGWAS
#' @return list of latent GWAS results in format for flashfmZero, with columns rsID, beta, EAF, beta_se, p-value
#' @export
#' @author Jenn Asimit
latentGWAS <- function(obsGWAS,covY,L,beta_colname="BETA",se_colname="SE",snpID_colname="SNP",EAF){
 
 L <- as.matrix(L)
 covY <- as.matrix(covY)
 Lterm <- LoadingTerm(L)
 P <- nrow(L) # number of observed traits
 K <- ncol(L) # number of latent factors
 
 beta_traits <- mapply(function(x) {x[,beta_colname]}, obsGWAS) # matrix of beta.hat, rows=snps, columns=traits
 se_traits <- mapply(function(x) {x[,se_colname]}, obsGWAS)
 
 betaL <- sapply(asplit(beta_traits, 1), latentBeta, Lterm=Lterm) # asplit makes a list where each component is a row of the matrix
 betaL <- t(betaL) # rows are variants, columns are traits

 seL <- sapply(asplit(se_traits, 1), latentBetaSE, Lterm = Lterm, covY = covY)
 seL <- t(seL)
 
 ZL <- betaL/seL
 
 pL <- 2*pnorm(-abs(ZL))
 
 colnames(betaL) <- colnames(seL) <- colnames(pL) <- colnames(L)
 rownames(betaL) <- rownames(seL) <- rownames(pL) <- obsGWAS[[1]][,snpID_colname]
 snpnames <- obsGWAS[[1]][,snpID_colname]
 
# return(list(betaL=betaL, seL=seL, pvalL=pL)) 
 gwasL.list <- vector("list",K)
 for(i in 1:K){
  gwasL.list[[i]] <- data.frame(rsID=snpnames, beta=betaL[,i], EAF=EAF, beta_se=seL[,i], p_value=pL[,i])
 }
 names(gwasL.list) <- colnames(L)

return(gwasL.list)
}






#' @title GWAS for latent factors at one variant using observed trait GWAS summary statistics
#' @param L PxK matrix of factor loadings generated from a factor analysis of P observed traits and K latent factors (K<P)
#' @return function of factor loadings (L^T L)^{-1}L^T; needed to estimate latent GWAS effect estimates from observed traits
#' @author Jenn Asimit
LoadingTerm <- function(L) {
 L <- as.matrix(L)
 mat <- t(L)%*%L
 matinv <- solve(mat)
 out <- matinv%*%t(L)
 return(out)
}



#' @title GWAS effect estimates for latent factors at one variant using observed trait GWAS summary statistics
#' @param beta_traits vector of effect estimates from a set of observed trait GWAS
#' @param Lterm loadingterm matrix output from LoadingTerm function
#' @return vector of effect estimates for the latent factors (at one variant)
#' @author Jenn Asimit
latentBeta <- function(beta_traits,Lterm) {
 betaL <- as.matrix(Lterm)%*%matrix(beta_traits,ncol=1)
 return(betaL)
}

#' @title GWAS effect estimate standard errors for latent factors at one variant using observed trait GWAS summary statistics
#' @param se_traits vector of standard errors 
#' @param Lterm loadingterm matrix output from LoadingTerm function
#' @param covY trait covariance matrix (same as used in factor analysis)
#' @return vector of effect estimate standard errors for the latent factors (at one variant)
#' @author Jenn Asimit
latentBetaSE <- function(se_traits,Lterm, covY) {
 var_beta_hat <- as.matrix(covY) * (matrix(se_traits, ncol=1) %*% matrix(se_traits, nrow=1)) # Var(beta1,beta2) = cov(Y1,Y2)*beta.se1*beta.se2
 varLa <- as.matrix(Lterm) %*% var_beta_hat
 varL <- varLa%*%t(Lterm)
 seL <- sqrt(diag(varL))
 return(seL)
}





#' @title Re-scale factor loadings to get their contributions and output factor loadings and contributions matrices with observed traits ordered by maximum contributing latent factor
#' @param L PxK matrix of factor loadings generated from a factor analysis of P observed traits and K latent factors (K<P)
#' @return factor loadings and re-scaled factor loadings (contributions) matrices with observed traits ordered by maximum contributing latent factor
#' @export
#' @author Jenn Asimit
factor_contributions <- function(L) {
 FAloading <- L
 prop.fn <- function(x) x^2/sum(x^2)
 PropLY <- t(apply(FAloading,1,prop.fn))
# maxLY <- unique(apply(PropLY,1,which.max)) # unique max contributing latent factors
# otherLY <- setdiff(1:ncol(PropLY),maxLY) # left-over latent factors
# colreorder <- c(maxLY,otherLY) 
# FAloading <- L[,colreorder]
# PropLY <- PropLY[,colreorder]
 obstrait_maxfactor <- apply(PropLY,1,which.max) # using new column re-ordering identify max contributing latent factor
 rowreorder <- names(sort(obstrait_maxfactor)) # rows appear by common max latent factors
 PropLY  <- PropLY[rowreorder,]
 FAloading <- FAloading[rowreorder,]
 return(list(rescaledFAloadings=PropLY, FAloading=FAloading)) 
}

