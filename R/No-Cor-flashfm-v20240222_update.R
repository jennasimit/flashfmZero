

#' @title Wrapper for flashfm Multi-Trait Fine-Mapping with JAM (when trait correlation is zero)  - this is the dynamic number of max causal variant version
#' @param gwas.list List of M data.frame objects, where M is the number of traits; gwas.list\[\[i\]\] is a data.frame for  trait i with 3 columns named: rsID, beta, EAF
#' @param corX SNP correlation matrix  
#' @param N Vector of length M; Nall\[i\] is the (effective) sample size for trait i
#' @param save.path Path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param TOdds target odds of no sharing to sharing; default is 1
#' @param cpp cumulative posterior probability threshold for selecting top models; default 0.99
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @param maxcv starting value for maximum number of causal variants; default 1
#' @param maxcv_stop maximum value to consider for maximum number of causal variants; maxcv_stop >= maxcv.
#' @param jam.nM.iter in millions, number of iterations to use in JAM; defailt 1 (1 million)
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.6
#' @param minsnpmppi only group snps with total MPPI > minsnpmppi; default 0.01
#' @return list with 2 components: mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups, 
#' a list with 2 components giving the SNP groups constructed under single-trait (snpGroups\[\[1\]\]) and multi-trait fine-mapping (snpGroups\[\[2\]\])
#' @author Jenn Asimit
# #' @import R2BGLiMS
#' @export
FLASHFMZEROwithJAMd <- function (gwas.list, corX, N, save.path, TOdds = 1, 
    cpp = 0.99, NCORES, maxcv = 1, maxcv_stop = 20, jam.nM.iter = 1, r2.minmerge=0.6, minsnpmppi = 0.01) 
{
    M <- length(gwas.list)
    ybar <- numeric(M)
 #   covY <- diag(1,M)
 #   Vy <- diag(covY)
 	Vy <- rep(1,M)
    corX <- as.matrix(corX)
    if (!dir.exists(save.path)) {
        message(c("Directory ", save.path, " does not exist. Creating directory ", 
            save.path))
        dir.create(save.path)
    }
    tmpdir <- paste0(save.path, "/tmp", sample(1:1000, 1))
    dir.create(tmpdir)
    main.input <- JAMmulti2_sameN(gwas.list, corX, ybar, Vy, N, r2 = 0.99, 
        save.path, maxcv = maxcv, maxcv_stop = maxcv_stop,  
        jam.nM.iter = jam.nM.iter,NCORES=NCORES)
    gc(verbose = FALSE)
#    ss.stats <- flashfm:::summaryStats(Xmat = FALSE, ybar.all = ybar, main.input = main.input)
	fm.multi <- flashfmZero(main.input, TOdds = TOdds, 
            cpp = cpp, maxmod = NULL, NCORES = NCORES)    
    snpGroups <- makeSNPgroups2U(main.input, fm.multi, is.snpmat = FALSE, 
        min.mppi = 0.01, minsnpmppi = minsnpmppi, r2.minmerge = r2.minmerge)
    #mpp.pp <- flashfm:::PPsummarise(fm.multi, snpGroups, minPP = 0.01)
    mpp.pp <- PPsummarise(fm.multi, snpGroups, minPP = 0.01)
    unlink(paste0(tmpdir, "/*"))
    return(list(mpp.pp = mpp.pp, snpGroups = snpGroups))
}


##' @title Marginal PP for models sharing information between traits, when trait correlation is zero 
##' @param STR list of models for traits 1, 2, ..., n, each given in
##'     the form of a character vector, with entries
##'     \code{"snp1\%snp2\%snp3"}. The null model is given by
##'     \code{"1"} OR \code{"0"}.  It is assumed that all elements of
##'     ABF, PP and pr below follow this same order.
##' @param PP list of posterir probabilities for the models in M
##' @param mbeta list of joint beta estimates for each trait
##' @param kappa single value or vector of values to consider for the
##'     sharing scale parameter.  the value of kappa=1 must be
##'     included, and if not will be prepended.
##' @param N number of individiduals with measurements for all traits
##' @param nsnps number of snps in region
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
##' @return list of: - single.pp: list of pp for each model in
##'     \code{STR[[i]]} for trait i - shared.pp: list of pp for each model
##'     in \code{STR[[i]]} for trait i
						
marginalpp0 <- function(STR, PP, mbeta, kappa, N,nsnps,NCORES) {  
    
    n <- length(STR) # number of traits
  
    if(n<2)
        stop("Need at least 2 traits")
    if( length(STR)!=n || length(PP)!=n )
        stop("STR and PP need to have the same lengths")
   
    if(is.null(names(STR)))
        names(STR) <- paste0("trait",seq_along(STR))
   qt <- names(STR)
    
    ## calculate model sizes 
    SS <- lapply(STR,strsplit,"%")
    usnps <- sort(unique(unlist(SS)))
    nsnpspermodel <- lapply(SS,function(x) sapply(x,length))
    for(i in seq_along(STR)) {
        wh <- which(STR[[i]] %in% c("0","1"))
        nsnpspermodel[[i]][wh] <- 0
    }
    maxsnps <- max(unlist(nsnpspermodel))
    tau <- outer(0:maxsnps,0:maxsnps,calctau,nsnps=nsnps,kappa=kappa)
   
  
   
    alt.pp <- calcAdjPP(qt=qt,STR=STR,SS=SS,tau=tau,nsnpspermodel=nsnpspermodel,kappa=kappa,PP=PP,beta=mbeta,NCORES)

    
    for(i in seq_along(alt.pp)){
 	names(alt.pp[[i]]) <- STR[[i]]
 	}
    ret <- lapply(seq_along(qt), function(i) {
        data.frame(single.pp=PP[[i]],
                   shared.pp=alt.pp[[i]])})
    names(ret) <- qt
    return(ret)
}


#' @title Marginal PP for models of a set of traits, sharing information between the traits, when trait correlation is zero
#' @param main.input List of 3 components: SM=list of snpmod objects for a set of traits; mbeta=list of joint effects for each trait; nsnps= number of SNPs in the region 
#' This could be obtained from flashfm.input or JAMexpanded.multi.
#' @param TOdds Vector of target odds of no sharing to sharing
#' @param cpp cumulative posterior probability threshold for selecting top models; this is ignored when maxmod is spespecified
#' @param maxmod maximum number of top models to output; NULL by default
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @return List consisting of PP: marginal PP for models and MPP: marginal PP of SNP inclusion
#' @export
#' @author Jenn Asimit
flashfmZero <- function(main.input,TOdds=1,cpp=0.99,maxmod=NULL,NCORES) {
	
	N <- main.input$Nall
	nsnps <- main.input$nsnps
	mbeta <- main.input$mbeta
	SM <- main.input$SM
	
	nd <- M <- length(SM)
	
	qt <- names(main.input$SM)    	
	kappas <- c()
	#for(j in 1:length(TOdds)) kappas <- c(kappas,flashfm:::calckappa(nsnps=nsnps,p=2/nsnps,ndis=nd,target.odds=TOdds[j]))
	for(j in 1:length(TOdds)) kappas <- c(kappas, calckappa(nsnps=nsnps,p=2/nsnps,ndis=nd,target.odds=TOdds[j]))
    kappas <- round(kappas)
    traits <- paste(qt, collapse = "-")
    bestmod.thr <- vector("list",M)
 for(i in 1:M) {
 	bm <- best.models.cpp(SM[[i]],cpp.thr=cpp,maxmod)   
 	bestmod.thr[[i]] <- bm$models
 	message("Trait ",i, " (",qt[i],") ", "has cpp before adjustment: ",bm$old.cpp)
 	}
   
 STR <- lapply(bestmod.thr, "[[", "str") 
 PP <- lapply(bestmod.thr, "[[", "PP")

 names(STR) <- qt
 names(PP) <- qt


 for(i in 1:M) mbeta[[i]] <- mbeta[[i]][STR[[i]]]

    pp <- vector("list",length=nd) 
       
     for(kappa in kappas) {
     
     ret <- marginalpp0(STR, PP, mbeta, kappa, N,nsnps, NCORES)    
     for(i in 1:nd) pp[[i]] <- cbind(pp[[i]],ret[[i]]$shared.pp)
     } 
      for(i in 1:nd) {
       pp[[i]] <- cbind(ret[[i]]$single.pp,pp[[i]])
       colnames(pp[[i]]) <- paste("pp",c("null",round(TOdds,2)),sep=".")
       rownames(pp[[i]]) <- rownames(ret[[i]])
       }

   
    mpp <- lapply(pp, MPP.fn)
    names(pp) <- qt
    mpp1 <- lapply(mpp, t)
   
    MPP <- mpp1[[1]] 
    for (k in 2:M) MPP <- gtools::smartbind(MPP, mpp1[[k]], fill = 0)
    return(list(PP = pp, MPP = MPP,sharing=c("null",kappas)))
}



#' @title Key input for flashfm - constructs snpmod object list and joint effect estimates list for all trait if have external single trait fine-mapping results, when trait correlation is zero
#' @param modPP.list list of data.frame objects for each trait, containing a column named "snps": snp models of the form \code{"snp1,snp2"} and "PP": snp model posterior probability from single trait fine-mapping
#' @param beta1.list list of single SNP effect estimates for each trait in the form of a named vector; name of each effect estimate should appear in Gmat
#' @param corX genotype matrix with SNP columns; could be from sample or reference panel; use unrelated samples
#' @param Nall vector of sample sizes; if related samples then give effective sample sizes
#' @param ybar.all vector of trait means; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set ybar as 0-vector
#' @param raf named vector of SNP reference allele frequencies where name is snp id (must match SNP coding used for effect estimates); only needed if Gmat is a covariance matrix and MUST be in same order as SNPs in covariance matrix
# #' @param is.snpmat is snpmat
#' @return list containing the main input for flashfm
#' @author Jenn Asimit
#' @export
#flashfm0.input <- function(modPP.list,beta1.list,corX,Nall,ybar.all,raf) {
flashfmZero.input <- function(modPP.list,beta1.list,corX,Nall,ybar.all, raf) {  
   is.snpmat=FALSE
   M <- length(Nall)
   Gmat <- corX
   if(is.null(rownames(Gmat)) | is.null(colnames(Gmat))) {
     if(is.null(raf)) stop("Both a covariance matrix and vector of SNP effect allele frequencies (raf) must be provided.")
     rownames(Gmat) <- colnames(Gmat) <- names(raf)
     }
   
   # filter snps to include only snps present for all traits and in snp reference/covariance matrix
	snps <- NULL
 	for(i in 1:M) snps <- union(snps,names(beta1.list[[i]]))
 	snps <- intersect(snps,colnames(Gmat))
 	for(i in 1:M) beta1.list[[i]] <- beta1.list[[i]][snps]
 	Gmat <- Gmat[snps,snps]; raf <- raf[snps] 
   
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
   

   
   return(list(SM=SM,mbeta=mbeta,nsnps=nsnps,Nall=Nall,Gmat=Gmat,beta1.list=beta1.list,raf=raf))
}



#' @title internal function for calcAdjPP for a pair of traits
#' @param i model index for trait 1
#' @param j model index for trait 2
#' @param T1 index of trait 1
#' @param T2 index of trait 2
#' @param SS list consisting of lists of model configuration SNPs for each trait
#' @param PP list consisting of lists of model PP for each trait
#' @param tau matrix of adjustment terms
#' @param nsnpspermodel list of number of SNPs per model for each model in STR
#' @param kappa single value of sharing parameter kappa
#' @author Jenn Asimit
calcQ12 <- function(i,j,T1,T2,SS,PP,tau,nsnpspermodel,kappa) {
# contributes to Q for 1 | 2 and 2|1
if(SS[[T1]][[i]][1] =="1" | SS[[T2]][[j]][1] == "1") { #at least one is null model -> tau=1 & intersection is empty
 kadj <- 1
 tij <- 1
 } else {
overlap <- 1*(any(SS[[T1]][[i]] %in% SS[[T2]][[j]]))  
kadj <- ifelse(overlap==0,1,kappa) 
tij <- tau[(nsnpspermodel[[T1]][i]+1),(nsnpspermodel[[T2]][j]+1)] # shift array indices by 1 since for numsnps 0 to maxnum
}
tk <- tij*kadj
adj1 <- PP[[T2]][[j]] * tk
adj2 <- PP[[T1]][[i]] * tk
return(c(adj1,adj2))
}

vcalcQ12 <- Vectorize(calcQ12,vectorize.args=c("i","j"),SIMPLIFY=FALSE) #last arg is so that have single element output and can apply outer






#' @title Calculates trait-adjusted posterior probabilities for all traits at sharing parameter kappa
#' @param qt vector of trait names
#' @param STR list consisting of vectors of model configurations for each trait
#' @param SS list consisting of lists of model configuration SNPs for each trait
#' @param tau matrix of adjustment terms
#' @param nsnpspermodel list of number of SNPs per model for each model in STR
#' @param kappa single value of sharing parameter kappa
#' @param PP list consisting of vectors of posterior probabilities for the model configurations for each trait
#' @param beta list of joint effect estimates for models in STR; multi.beta output
# #' @param SSy matrix of trait cross-products
# #' @param Sxy matrix with each column being the cross-product between SNPs and a trait
# #' @param xcovo SNP covariance matrix
# #' @param Mx vector of SNP means
# #' @param N number of individuals with measurements for all traits
# #' @param allVres list of variance residuals
# #' @param covY covariance matrix of traits
# #' @param Nqq matrix of all pair-wise counts of number of individuals with both traits in a pair measured;
# #' @param Nq3  vector of counts of number of individuals with three traits measured; all triples considered; NULL if M < 4
# #' @param Nq4  vector of counts of number of individuals with four traits measured; all quadruples considered; NULL if M < 5
# #' @param fastapprox logical that is TRUE when fast approximation is used that does not include unequal sample size adjustments; default is FALSE
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1; 
#' @return list of trait-adjusted posterior probabilities for each trait at sharing parameter kappa
#' @author Jenn Asimit
calcAdjPP <- function(qt,STR,SS,tau,nsnpspermodel,kappa,PP,beta,NCORES) {
 
    M <- length(qt)
    np <- choose(M,2)
    c2 <- combn(1:M,2,simplify=TRUE)
    c2names <- apply(c2,2, function(cc) return(paste0("Q",paste(cc,collapse=".Q"))))
	Q <- structure(vector("list",np),names=c2names)
	nummods <- sapply(STR,length)

	for(i in 1:np) { # for each qt pair Q[[i]] is a matrix where Q[[i]][j,k] is a list with 
					# two components adjPP12[modj for T1,modk for T2], adjPP21[modj for T1,modk for T2] where 1=c2[1,i], 2=c2[2,i]
	    
     tmp <- outer(1:nummods[c2[1,i]],1:nummods[c2[2,i]],vcalcQ12,T1=c2[1,i],T2=c2[2,i],SS,PP,tau,nsnpspermodel,kappa)
     q12 <- apply(tmp,2,function(x) unlist(lapply(x,"[[",1))) # first element of each cell in matrix; q12 is a matrix
 	 q21 <- apply(tmp,2,function(x) unlist(lapply(x,"[[",2)))  	    
	 q1 <- log(q12)
     q2 <- t(log(q21))
     q1 <- apply(q1,1,logsum); q1 <- q1-logsum(q1)
     q2 <- apply(q2,1,logsum); q2 <- q2-logsum(q2)
     Q[[i]] <- list(q1,q2)   # q1 and q2 are vectors, possibly different lengths
		}
	
	qns <- unlist(strsplit(names(Q),"[.]"))
	PPadj <- structure(vector("list",M),names=qt)
	
	if(M==2) { 
	  i=1#	  delta <- calcD12(1:nummods[c2[1,i]],1:nummods[c2[2,i]],T1=c2[1,i],T2=c2[2,i],beta=beta,SSy=SSy,Sxy=Sxy,xcovo=xcovo,Mx=Mx,Nqq=Nqq,Vres=allVres,covY=covY,nsnpspermodel)
      tmp <- Q[[i]]
      q12 <- Q[[i]][[1]]
 	  q21 <- Q[[i]][[2]]	    
 	  
 	  lPP <- lapply(PP,log)
 	  PPadj[[1]] <- lPP[[1]] + q12 
 	  PPadj[[2]] <- lPP[[2]] + q21
 	  PPadj[[1]] <- exp(PPadj[[1]] - logsum(PPadj[[1]]))
 	  PPadj[[2]] <- exp(PPadj[[2]] - logsum(PPadj[[2]]))
 	  names(PPadj) <- qt	
#	  PPadj[[1]] <- PP[[1]]*q1/sum(PP[[1]]*q1)
#	  PPadj[[2]] <- PP[[2]]*q2/sum(PP[[2]]*q2)
	  
	} else {
		lPP <- lapply(PP,log)
		
		PPadj <- vector("list",M)
		ivec <- vector("list",M)
		for(i in 1:M) ivec[[i]] <- i
		Tadj <- parallel::mclapply(ivec,pre.ppadj,qns,Q,mc.cores =NCORES)	
		for(i in 1:M) {
			PPadj[[i]] <- lPP[[i]] + Tadj[[i]] 	
			PPadj[[i]] <- exp(PPadj[[i]] - logsum(PPadj[[i]]))
			}
		names(PPadj) <- qt	
     	}
	return(PPadj) 
	 
}


#### ppadj functions ####

pre.ppadj <- function(i,qns,Q) {
	    
	    qn  <- paste0("Q",i)
	 	ind <- grep(qn,qns,fixed=TRUE)
	 	whO <- ind[which(ind %% 2 == 1)] # odd indices so first list component 	 
	 	whE <- ind[which(ind %% 2 == 0)]
	 	keep <- NULL
	 	if(length(whO)>0) {
	 	  tmp <- Q[(whO+1)/2]
	 	  keep <- lapply(tmp, function(x) x[[1]])
	 	  nk <- names(keep)
	 	  names(keep) <- unlist(strsplit(nk,"[.]"))[c(FALSE,TRUE)]
	 		}
	 	if(length(whE)>0) {
	 	 tmp <- Q[whE/2]
	 	 if(!is.null(keep)) {
	 	 keep2 <- lapply(tmp,function(x) x[[2]] )
	 	 nk <- names(keep2)
	 	 names(keep2) <- unlist(strsplit(nk,"[.]"))[c(TRUE,FALSE)]
	 	 keep <- append(keep,keep2)	# 2nd component in list pair
	 	 				} else {
	 	 				keep <- lapply(tmp,function(x) x[[2]] )
	 	 				nk <- names(keep)
	 	 				names(keep) <- unlist(strsplit(nk,"[.]"))[c(TRUE,FALSE)]
	 	 				}
	 	} 
	 	keep <- as.matrix(do.call(cbind, keep))
	 	out <- apply(keep,1,sum) # log-scale
	 	return(out)
	 	}
