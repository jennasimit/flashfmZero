
#' @title internal processing function for JAMexpanded.multi
#' @param binout named vector of indicators for SNP presence in model; names are snp ids
#' @return snp model represented by snps separated by \code{"\%"}
#' @author Jenn Asimit
mod.fn <- function(binout) {
 msnps <- names(binout)
 ind <- which(binout==1)
 if(length(ind) > 1) { 
 m <-  paste(msnps[ind],collapse="%")
 } else if(length(ind) == 1) { 
   m <- msnps[ind]
  } else {m <- "1"}
return(m)
}

#' @title internal function for expanding models by tag SNPs in JAMexpanded.multi
#' @param snpmod intial snp model with snps separated by \code{"\%"}
#' @param taglist list of tag snps for each snp
#' @author Jenn Asimit
expand.mod<- function(snpmod,taglist) {
 if(snpmod == "1") {
  df <- data.frame(str=snpmod,snps=snpmod,size=0,tag=FALSE,stringsAsFactors = FALSE)
 } else { 
 snps <- unlist(strsplit(snpmod,"%"))
 ns <- length(snps)
 tsnps <- vector("list",ns)
 for(i in 1:ns) {
  tmp <- taglist[taglist$SNP==snps[i],"TAGS"]
  if(tmp == "NONE") {
    tsnps[[i]] <- snps[i]
     } else {
      tsnps[[i]] <- c(snps[i],unlist(strsplit(tmp,"[|]"))) 
      }
  }    
 emods <- expand.grid(tsnps,stringsAsFactors = FALSE)  
 out <- apply(emods,1,function(x){paste0(x,collapse="%")})
 Imod <- which(out==snpmod)
 istag <- rep(TRUE,length(out))
 istag[Imod] <- FALSE
 df <- data.frame(str=snpmod,snps=out,size=ns,tag=istag,stringsAsFactors = FALSE)
 
 } 
 
 return(df)
}

#' @title internal function for multibeta, modified from JAM_PointEstimates_Package_Simplified of R2BGLiMS package (Paul Newcombe)
#' @param marginal.betas vector of single snp effect estimates from GWAS
#' @param X.ref genotype matrix 
#' @param n sample size of GWAS
#' @param ybar trait mean
#' @return joint effect estimates of snps
JAM_PointEstimates_updated <- function(marginal.betas=NULL,X.ref=NULL,n=NULL,ybar) {
  
  # --- Setup sample sizes
  n.ref <- nrow(X.ref)
  if (is.null(n)) {
    n <- n.ref
  }
  
  ################################################
  # --- Construct the z = X'y outcome vector --- #
  ################################################
  # The element for each SNP is constructed from:
  # 1) Infer predicted y-values from the marginal.betas. Then mean-center for 0-intercept model
  # 2) Matrix multiply X.ref by the predicted y-values
  z <- rep(NA,length(marginal.betas))
  names(z) <- names(marginal.betas)
  for (v in 1:length(marginal.betas)) {
    y.pred <- X.ref[,v]*marginal.betas[v]
    y.pred.centered <- y.pred - mean(y.pred)
    z[v] <- X.ref[,v] %*% y.pred.centered # t(X.ref)%*%y
  }
  z <- z*n/n.ref
  
  #################################################
  # --- Construct multivariate beta estimates --- #
  #################################################
  
  # 1) Mean-centre X.ref
  xbar <- apply(X.ref,2,mean)
  for (v in 1:ncol(X.ref)) {
    X.ref[,v] <- X.ref[,v] - xbar[v] # MUST mean-center since z is constructed under 0 intercept 
  }
  
  # 2) Calculate MLE corresponding to the summary model
  xtx <- (t(X.ref) %*% X.ref)*n/n.ref
  multivariate.beta.hat <- solve(xtx) %*% z
  
  b0 <- ybar-sum(xbar*multivariate.beta.hat)
  out <- rbind(b0,multivariate.beta.hat)
  rownames(out) <- c("one",names(marginal.betas))

  return(out)
}

#' @title internal function for multibeta, modified from JAM_PointEstimates_Package_Simplified of R2BGLiMS package
#' @param marginal.betas vector of single snp effect estimates from GWAS
#' @param Xcov genotype covariance matrix
#' @param raf vector of reference allele frequencies 
#' @param n sample size of GWAS
#' @param ybar trait mean
#' @return joint effect estimates of snps
#' @author Jenn Asimit
JAM_PointEstimates_Xcov <- function(marginal.betas=NULL,Xcov=NULL,raf,n,ybar) {
 
  xbar <- 2*raf
  
  
  ################################################
  # --- Construct the z = X'y outcome vector --- #
  ################################################
  # The element for each SNP is constructed from:
  # 1) cov(x,y) = beta.hat[x]*V(x)
  # 2) x'y = n*E[XY] = n*(cov(x,y)) under zero intercept model
  z <- rep(NA,length(marginal.betas))
  names(z) <- names(marginal.betas)
  for (v in 1:length(marginal.betas)) {
	z[v] <- n*(marginal.betas[v]*Xcov[v,v])
  }
  
  
  #################################################
  # --- Construct multivariate beta estimates --- #
  #################################################
  # under zero intercept model
  
  xtx <- n*(Xcov)
  multivariate.beta.hat <- as.vector(solve(xtx) %*% z)
  
  
  # add in intercept
  b0 <- ybar-sum(as.vector(xbar)*multivariate.beta.hat)
  out <- matrix(c(b0,multivariate.beta.hat),ncol=1)
  rownames(out) <- c("one",names(marginal.betas))

  return(out)
}


#' @title Using summary statistics, calculates joint effect estimates
#' @param mod joint SNP model with snps separated by \code{"\%"} e.g. \code{"snp1\%snp2"}
#' @param beta1 named vector single-SNP effect estimates; each effect estimate has name given by SNP id and ids must appear in Gmat columns
#' @param Gmat genotype matrix with SNP columns and individuals rows; could be reference matrix or from sample
#' @param N sample size for trait 
#' @param ybar trait mean
#' @param is.snpmat logical taking value TRUE when Gmat is a genotype matrix and FALSE when Gmat is a SNP covariance matrix
#' @param raf named vector of SNP reference allele frequencies (must match SNP coding used for effect estimates); only needed if Gmat is a covariance matrix
#' @return joint effect estimates for SNPs in model given by mod
#' @author Jenn Asimit
#' @export
multibeta <- function(mod,beta1,Gmat,N,ybar,is.snpmat,raf=NULL) {

 if(mod == "1") {
  mbeta <- ybar; names(mbeta) <- "one"
 } else {
  msnps <- unlist(strsplit(mod,"%"))
  if(all(msnps %in% names(beta1))) { 
   B <- beta1[msnps]
   if(is.snpmat) {  
     if(all(msnps %in% colnames(Gmat))){ mbeta <- JAM_PointEstimates_updated(B,X.ref=as.matrix(Gmat[,msnps]),n=N,ybar=ybar)
       } else { mbeta <- NA}
     } else { 
 	       colnames(Gmat) <- names(raf)
 	       rownames(Gmat) <- names(raf)
 	       if(all(msnps %in% colnames(Gmat))) { mbeta <- JAM_PointEstimates_Xcov(B,Xcov=as.matrix(Gmat[msnps,msnps]),raf=raf[msnps],n=N,ybar=ybar)  
 	       } else { mbeta <- NA}
 	       }
   } else {mbeta <- NA }
 }
 
 return(mbeta)
} 

 
#' @title Calculate approximate Bayes' factor (ABF) 
#' @param mod joint SNP model with snps separated by \code{"\%"} e.g. \code{"snp1\%snp2"}
#' @param mbeta joint effect estimates for SNPs in model given by mod; output from multibeta
#' @param SSy sum(y.squared) for trait y
#' @param Sxy vector of sum(xy) for snp x, trait y
#' @param Vy trait variance
#' @param N sample size
#' @return ABF for model mod
#' @author Jenn Asimit
#' @export
calcABF <- function(mod,mbeta,SSy,Sxy,Vy,N) {
 if(mod=="1") { msnps <- "one"
 } else {msnps <- c("one",unlist(strsplit(mod,"%")))}
 beta <- mbeta[[mod]]
 num <- SSy-sum(Sxy[msnps]*beta)
 den <- (N-1)*Vy
 k <- length(msnps)-1
 out <- -N*.5*log(num/den)-0.5*k*log(N)
 return(out)
}

# SNP prior based on binomial distribution; modification of snpprior function (Chris Wallace)
SNPprior <- function(x=0:10, n, expected, overdispersion=1, pi0=NA, truncate=NA, overdispersion.warning=TRUE) {
  if(overdispersion < 1 & overdispersion.warning)
    stop("overdispersion parameter should be >= 1")
  x <- as.integer(x)
  if(any(x<0))
    stop("x should be an integer vector >= 0")
  if(!is.na(truncate))
    x <- seq(min(x), max(min(truncate,n),x))
  if(any(x>n))
    stop("max x should be <= n")
  p <- expected/n
    rho <- (overdispersion - 1)/(n-1)
    
	prob <- log(dbinom(x, size=n, prob=p))
#    prob <- if(rho==0) {
#                log(dbinom(x, size=n, prob=p))
#            } else {
#                log(dbetabinom(x, size=n, prob=p, rho=rho))
#            }

  if(!is.na(pi0) && 0 %in% x)
    prob[x==0] <- log(pi0)
  if(!is.na(truncate))
    prob <- prob - logsum(prob)
  prob <- prob - lchoose(n,x)  
  names(prob) <- as.character(x)
  return(exp(prob))
}

# converts data.frame of abfs to snpmod object; modified abf2snpmod (Chris Wallace)
makesnpmod <- function (abf, expected, overdispersion = 1, nsnps = NULL) {
    tmp <- new("snpmod")
    msize <- nchar(gsub("[^%]", "", abf$model)) + 1
    msize[abf$model == "1"] <- 0
    
    prior <- SNPprior(x = 0:max(msize), expected = expected, 
        n = nsnps, truncate = max(c(msize, 20)), overdispersion = overdispersion)
    mprior <- prior[as.character(msize)]
    mpp <- log(mprior) + abf$lBF
    mpp <- mpp - logsum(mpp)
    message("creating snpmod data.frame")
    tmp@models <- data.frame(str = abf$model, logABF = abf$lBF, 
        by.tagging = abf$tag, size = msize, prior = mprior, 
        lPP = mpp, PP = exp(mpp), stringsAsFactors = FALSE)
    tmp@model.snps <- strsplit(tmp@models$str, "%")
    message("calculating marginal SNP inclusion probabilities")
    marg.snps(tmp)
}

#' internal function marg.snps 
#' @param d snpmod object
#' @author Chris Wallace
marg.snps <- function(d) {
    mod <- makemod(d@models$str)
    marg.pp <- (d@models$PP %*% mod)[1,,drop=TRUE]
    d@snps <- marg.snps.vecs(d@models$str, d@models$PP)
    return(d)
}

#' internal function makemod 
#' @param snps model given by snp ids separated by \code{"\%"}
#' @author Chris Wallace
makemod <- function(snps) {
  snps <- strsplit(snps,"%")
  all.snps <- unique(unlist(snps))
  mod <- Matrix::Matrix(0,length(snps),length(all.snps),dimnames=list(NULL,all.snps),doDiag=FALSE)
  snum <- lapply(snps, function(s) which(all.snps %in% s))
  I <- rep(1:length(snum),times=sapply(snum,length))
  J <- unlist(snum)
  mod[cbind(I,J)] <- 1
  mod
}

#' internal function marg.snps.vec 
#' @param str models given by snp ids separated by \code{"\%"}
#' @param pp posterior probabilities
#' @author Chris Wallace
marg.snps.vecs <- function(str,pp) {
    mod <- makemod(str)
    marg.pp <- (pp %*% mod)[1,,drop=TRUE]
    data.frame(Marg_Prob_Incl=marg.pp, var=names(marg.pp), rownames=names(marg.pp), stringsAsFactors=FALSE)
}


#' @title Key input for flashfm - constructs snpmod object list and joint effect estimates list for all trait if have external single trait fine-mapping results
#' @param modPP.list list of data.frame objects for each trait, containing a column named "snps": snp models of the form \code{"snp1,snp2"} and "PP": snp model posterior probabliity from single trait fine-mapping
#' @param beta1.list list of single SNP effect estimates for each trait in the form of a named vector; name of each effect estimate should appear in Gmat
#' @param Gmat genotype matrix with SNP columns; could be from sample or reference panel; use unrelated samples
#' @param Nall vector of sample sizes; if related samples then give effective sample sizes
#' @param ybar.all vector of trait means; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set ybar as 0-vector
#' @param related logical indicating if samples are related (TRUE) or not (FALSE); default is FALSE
#' @param y (optional) matrix of trait values (trait columns) or indicators of trait measured; used to get joint sample counts; default is NULL and if not provided an approximation is used based on vector of trait sample sizes
#' @param Nsame (optional) single sample size value, if all traits measured on all individuals
#' @param is.snpmat logical taking value TRUE when Gmat is a genotype matrix and FALSE when Gmat is a SNP covariance matrix
#' @param raf named vector of SNP reference allele frequencies where name is snp id (must match SNP coding used for effect estimates); only needed if Gmat is a covariance matrix and MUST be in same order ans snps in covariance matrix
#' @return list containing the main input for flashfm
#' @author Jenn Asimit
#' @export
flashfm.input <- function(modPP.list,beta1.list,Gmat,Nall,ybar.all,related=FALSE,y=NULL,Nsame=NULL,is.snpmat,raf=NULL) {
   M <- length(Nall)
   if(!is.snpmat) {
     if(is.null(raf)) stop("When is.snpmat=FALSE, both a covariance matrix and vector of SNP RAFs must be provided.")
     rownames(Gmat) <- colnames(Gmat) <- names(raf)
     }
   
   # filter snps to include only snps present for all traits and in snp reference/covariance matrix
	snps <- NULL
 	for(i in 1:M) snps <- union(snps,names(beta1.list[[i]]))
 	snps <- intersect(snps,colnames(Gmat))
 	for(i in 1:M) beta1.list[[i]] <- beta1.list[[i]][snps]
 	if(is.snpmat) {Gmat <- Gmat[,snps]
 	} else{ Gmat <- Gmat[snps,snps]; raf <- raf[snps] }
   
   mbeta <- vector("list",M)
   SM <- vector("list",M) 
   
   for(i in 1:M) {
    modsnps <- strsplit(modPP.list[[i]]$snps, ",")
    snpmods <- sapply(modsnps,function(x) paste(x,collapse="%"))
    modPP.list[[i]]$snps <- snpmods
    modPP.list[[i]]$str <- snpmods
    mbeta[[i]] <- lapply(modPP.list[[i]]$snps,multibeta,beta1.list[[i]],Gmat,Nall[i],ybar.all[i],is.snpmat=is.snpmat,raf=raf)	
    names(mbeta[[i]]) <- modPP.list[[i]]$snps
    mod.keep <- which(!is.na(mbeta[[i]]))
    mbeta[[i]] <- mbeta[[i]][mod.keep]
    ppdf <- modPP.list[[i]][mod.keep,c("str","PP")]
    ppdf$PP <- ppdf$PP/sum(ppdf$PP) # adjust to account for any removed models, so that PPs still sum to 1 
    SM[[i]] <- PP2snpmod(ppdf)
   } 
   nsnps <- ncol(Gmat)
   names(SM) <- names(modPP.list)
   names(mbeta) <- names(modPP.list)
   
   	if(!related) Nlist <- makeNlist(Nall,y,Nsame)
	if(related) Nlist <- makeNlist.rel(Nall,y,Nsame)
   
   if(is.null(raf)) raf <- apply(Gmat,2,mean)
   
   return(list(SM=SM,mbeta=mbeta,nsnps=nsnps,Nlist=Nlist,Gmat=Gmat,beta1.list=beta1.list,raf=raf))
}

#' @title Summary statistics needed for flashfm input
#' @param Xmat logical, TRUE if main.input is based on genotype matrix (reference or from sample); FALSE if main.input based on covariance matrix 
#' @param ybar.all vector of trait means
#' @param main.input output from flashfm.input
#' @return list of 4 components: Mx = mean of SNPs, xcovo = covariance matrix of SNPs including column and row of 0s for intercept term, Sxy = matrix of Sxy values (column traits), ybar=vector of trait means
#' @author Jenn Asimit
#' @export
summaryStats <- function(Xmat=TRUE,ybar.all,main.input) {

 Nall <- diag(main.input$N$Nqq)
 Xinfo <- main.input$Gmat
 Xmean <- 2*main.input$raf
 beta1.list <- main.input$beta1.list

  if(Xmat==TRUE) {
  Mx <- c("one"=1,apply(Xinfo,2,mean))
  xcov <- var(Xinfo)
  xcovo <- matrix(0,nrow=length(Mx),ncol=length(Mx),dimnames=list(names(Mx),names(Mx)))
  xcovo[2:ncol(xcovo),2:ncol(xcovo)] <- xcov
  xcovo[,1] <- 0; xcovo[1,] <- 0
  Vx <- diag(xcov)
  } else {
  Mx <- c("one"=1, Xmean)
  xcovo <- matrix(0,nrow=length(Mx),ncol=length(Mx),dimnames=list(names(Mx),names(Mx)))
  xcovo[2:length(Mx), 2:length(Mx)] <- as.matrix(Xinfo)
  xcovo[,1] <- 0; xcovo[1,] <- 0
  Vx <- diag(as.matrix(Xinfo))
  }
  
  M <- length(beta1.list)
  Sxy <- NULL
  for(i in 1:M) Sxy <- cbind(Sxy, c(Sxy.hat(beta1=beta1.list[[i]],Mx=Mx[-1],N=Nall[i],Vx=Vx,muY=ybar.all[i]),"1"=ybar.all[i]*Nall[i]))
  rownames(Sxy)[nrow(Sxy)] <- "one"
  colnames(Sxy) <- names(beta1.list)
  return(list(Mx=Mx, xcovo=xcovo, Sxy=Sxy,ybar=ybar.all))
}


# internal function for groupIDs.fn
ingroup.fn <- function(snp,group) {
 ind <- grep(snp,group)
 1*(length(ind)>0)
 }


#' @title Find SNP group ids for a set of SNPs
#' @param snpgroups  list of snpgroups
#' @param Msnps vector of snps
#' @return vector same length of msnps that gives the SNP group containing each SNP, or the SNP id if it does not belong to a group
#' @author Jenn Asimit
#' @export
groupIDs.fn <- function(snpgroups,Msnps) {
# outputs snp group for each snp in Msnps; if not in a group, output rsid
 ng <- length(snpgroups)
 ns <- length(Msnps)
 
 G <- character(ns)
 for(i in 1:ns) {
  if(Msnps[i]=="1") {
    G[i] <- "null"
  } else {
  check <- unlist(lapply(snpgroups,ingroup.fn,snp=Msnps[i]))
  if(sum(check)>0) { G[i] <- names(which(unlist(check)>0))
  } else {
          	G[i] <- Msnps[i]
                }
   }
  }
 return(G)
}


