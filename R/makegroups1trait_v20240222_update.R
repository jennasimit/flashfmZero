
#' @title Wrapper for JAMdynamic single-trait fine-mapping applied to a set of traits, and output SNP groups
#' @param gwas.list List of M data.frames, where M is the number of traits; gwas.list\[\[i\]\] is a data.frame (one for each trait) 
#' with 3 columns named: rsID, beta, EAF; 
#' if trait names are provided for the M data.frames, these trait names are given in output
#' @param corX SNP correlation matrix that must have snp names in row.names and col.names
#' @param N scalar that is the trait sample size - same for all traits; if different sample sizes across traits, use flashfm to return single and multi-trait results
#' @param save.path Path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param maxcv starting value for maximum number of causal variants
#' @param maxcv_stop maximum value to consider for maximum number of causal variants; maxcv_stop >= maxcv.
#' @param jam.nM.iter in millions, number of iterations to use in JAM; default 1 (1 million)
#' @param min.mppi trim snp groups with total MPPI < min.mppi in all diseases; default 0.01
#' @param minsnpmppi only group snps with total MPPI > minsnpmppi; default 0.01
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.6
#' @param NCORES number of cores for parallel computing; recommend NCORES=M if running on hpc, but if on Windows/Mac, use NCORES=1
#' @param extra.java.arguments A character string to be passed through to the java command line. E.g. to specify a
#' different temporary directory by passing "-Djava.io.tmpdir=/Temp".
#' @return List consisting of two objects: mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups, 
#' a list giving the SNP groups construced under single-trait fine-mapping
#' @export
#' @author Jenn Asimit
multiJAMd <- function (gwas.list, corX, N, save.path, 
     maxcv = 1, maxcv_stop = 20, 
    jam.nM.iter = 1, min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.6,NCORES=1,extra.java.arguments=NULL) 
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
        jam.nM.iter = jam.nM.iter,NCORES=NCORES,extra.java.arguments=extra.java.arguments)
    gc(verbose = FALSE)
#    ss.stats <- flashfm:::summaryStats(Xmat = FALSE, ybar.all = ybar, main.input = main.input)
#	fm.multi <- flashfm0(main.input, TOdds = TOdds, 
#            cpp = cpp, maxmod = NULL, NCORES = NCORES)    
    snpGroups <- makeSNPgroupsU(main.input, is.snpmat = FALSE, 
        min.mppi = min.mppi, minsnpmppi=minsnpmppi, r2.minmerge = r2.minmerge)
    mpp.pp <- lapply(main.input$SM, function(x) {ppdf <- x@models[,"PP",drop=F]; rownames(ppdf)<- x@models$str; PPsummariseST(ppdf,snpGroups)})    
    out <- list(MPP=lapply(mpp.pp,"[[","MPP"), MPPg=lapply(mpp.pp,"[[","MPPg"), PP=lapply(mpp.pp,"[[","PP") , PPg=lapply(mpp.pp,"[[","PPg"))
#    mpp.pp <- flashfm:::PPsummarise(fm.multi, snpGroups, minPP = 0.01)
	


    unlink(paste0(tmpdir, "/*"))
    return(list(mpp.pp =out, snpGroups = snpGroups))
}


#' @title Make SNP groups using fine-mapping information from all of the traits
#' @param main.input output from flashfm.input function
#' @param is.snpmat logical taking value TRUE when genotype matrix is provided and FALSE when covariance matrix is given
#' @param min.mppi trim snp groups with total MPPI < min.mppi in all diseases; default 0.01
#' @param minsnpmppi only group snps with total MPPI > minsnpmppi; default 0.01
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.5
#' @return list of  SNP groups
#' @export
makeSNPgroupsU <- function (main.input, is.snpmat, min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.5) 
{
    snp.data <- main.input$Gmat
    SMlist <- main.input$SM
    Xmat <- as.matrix(snp.data)
    sg <- groupmultiU(SMlist, Xmat, is.snpmat, min.mppi, minsnpmppi = minsnpmppi, 
        r2.minmerge)
    if (length(sg) > 1) {
        snpgroups <- sg$groups@.Data
        ng <- length(snpgroups)
        names(snpgroups) <- LETTERS[1:ng]
        if (ng > 26) 
            names(snpgroups)[27:min(ng, 52)] <- paste0(LETTERS[1:min(26, 
                ng - 26)], 2)
        if (ng > 52) {
            names(snpgroups)[27:min(ng, 52)] <- paste0(LETTERS[1:min(26, 
                ng - 26)], 2)
            names(snpgroups)[53:min(ng, 78)] <- paste0(LETTERS[1:min(26, 
                ng - 52)], 3)
        }
        if (ng > 78) {
            names(snpgroups)[27:min(ng, 52)] <- paste0(LETTERS[1:min(26, 
                ng - 26)], 2)
            names(snpgroups)[53:min(ng, 78)] <- paste0(LETTERS[1:min(26, 
                ng - 52)], 3)
            names(snpgroups)[79:min(ng, 104)] <- paste0(LETTERS[1:min(26, 
                ng - 78)], 4)
        }
        Ng <- lapply(snpgroups, length)
        sgd <- t(data.frame(Ng))
        colnames(sgd) <- "Group Size"
        message("SNP group sizes are: ")
        print(sgd)
    }
    else {
        snpgroups <- sg
    }
    return(snpgroups)
}



#' @title Return credible sets from single-trait fine-mapping of several traits
#' @param multiJAMout output from the function multiJAMd
#' @param cred probability for credible set; default is 0.99.
#'
#' @return List of two items: fm is a list of single-trait credible sets for several traits; cred is the probability for constructing the credible sets. 
#' @export
#' @author Jenn Asimit
multiJAMdCS <- function (multiJAMout, cred = 0.99) {
    mpp.pp <- multiJAMout$mpp.pp
    M <- length(mpp.pp$PP)
    csfm <-  vector("list", M)
    for (i in 1:M) {
    	ppmod <- mpp.pp$PP[[i]][,"PP"]
    	names(ppmod) <- rownames(mpp.pp$PP[[i]])
        #cs <- MGflashfm:::credsetU(ppmod, cred)
        cs <- credsetU(ppmod, cred)
        mpp <- mpp.pp$MPP[[i]][cs, ]
        csfm[[i]] <- data.frame(SNP = cs, MPP = mpp)
        csfm[[i]] <- csfm[[i]][order(csfm[[i]]$MPP, decreasing = TRUE), ]
    }
    names(csfm) <- names(mpp.pp$PP)
    return(list(fm = csfm, cred = cred))
}



JAMmulti2_sameN <- function(gwas.list, corX, ybar, Vy, N, r2 = 0.99, save.path, 
    maxcv = 2, maxcv_stop = 20, jam.nM.iter = 1, NCORES=1,extra.java.arguments=NULL) 
{
    maxcv_autocheck = TRUE
    M <- length(gwas.list)
    if (is.null(names(gwas.list))) {
        ts <- paste0("T", 1:M)
        names(ybar) <- ts
    }
    else {
        ts <- names(gwas.list)
    }
    beta1 <- lapply(gwas.list, function(x) {
        b <- x[, "beta"]
        names(b) <- x[, "rsID"]
        b
    })
    names(beta1) <- ts
    raf1 <- lapply(gwas.list, function(x) {
        b <- x[, "EAF"]
        names(b) <- x[, "rsID"]
        b
    })
    names(raf1) <- ts
#    Nlist <- makeNlist.rel(Ne = N)
#    N <- diag(Nlist$Nqq)
    snps <- Reduce(intersect, lapply(beta1, names))
    snps <- intersect(snps, colnames(corX))
    for (i in 1:M) {
        beta1[[i]] <- beta1[[i]][snps]
        raf1[[i]] <- raf1[[i]][snps]
    }
    corX <- as.matrix(corX)
    corX <- corX[snps, snps]
    nsnps <- length(snps)
    corX <- as.matrix(corX)
    #reftags <- MGflashfm:::cor.refdata2(corX, r2)
    reftags <- cor.refdata2(corX, r2)
    refGt <- as.matrix(reftags$refG)
    taglist <- reftags$taglist
    refG <- lqmm::make.positive.definite(refGt)
    refG <- as.matrix(refG)
    out <- list(SM = NULL, mbeta = NULL, Nlist = N)
    out$SM <- vector("list", M)
    names(out$SM) <- ts
    out$mbeta <- vector("list", M)
    names(out$mbeta) <- ts
    dd <- vector("list", M)
    SSy <- vector("list", M)
    Sxy <- vector("list", M)
    
#    for (j in 1:M) { outJ[[j]] <-  Jamext(j, raf=raf1[[j]], beta1=beta1[[j]], corX, Vy=Vy[j], refG, save.path, maxcv, maxcv_stop, jam.nM.iter, ybar=ybar[j], N )
#		}
 
 	ivec <- vector("list", M)
    for (i in 1:M) ivec[[i]] <- i
    outJ <- parallel::mclapply(ivec, JAMext, raf1, beta1, corX, refG, save.path, maxcv, maxcv_stop,jam.nM.iter, N, taglist,nsnps, mc.cores = NCORES,extra.java.arguments=extra.java.arguments)
    names(outJ) <- ts
    
    mbeta <- SM <- vector("list",M)
    for(i in 1:M) {
    	mbeta[[i]] <- outJ[[i]]$mbeta
    	SM[[i]] <- outJ[[i]]$SM
    }
    
    out$mbeta <- mbeta
    out$SM <- SM       
    out$N <- N
    out$nsnps <- nsnps
    #out$Gmat = flashfm:::cor2cov(corX, sd = sqrt(2 * raf1[[1]] * (1 - raf1[[1]])))
    out$Gmat = cor2cov(corX, sd = sqrt(2 * raf1[[1]] * (1 - raf1[[1]])))
    out$beta1.list = beta1
    out$raf = raf1[[1]]
    return(out)
}


JAMext <- function(j, raf1,beta1,corX, refG, save.path, maxcv, maxcv_stop, jam.nM.iter, N, taglist,nsnps,extra.java.arguments=NULL ){
 		ybar=0
 		Vy=1
 		raf <- raf1[[j]]
        BETA <- beta1[[j]][colnames(refG)]
        maf <- raf * (raf <= 0.5) + (1 - raf) * (raf > 0.5)
        mafs.ref <- maf[colnames(refG)]
        #covX <- flashfm:::cor2cov(corX, sd = sqrt(2 * raf * (1 - raf)))
        covX <- cor2cov(corX, sd = sqrt(2 * raf * (1 - raf)))
        covX <- as.matrix(covX)
        #topmods <- MGflashfm:::JAM.tries.maxcv(BETA, Vy, refG, mafs.ref, 
        topmods <- JAM.tries.maxcv(BETA, Vy, refG, mafs.ref,                                        
            N, save.path, maxcv = maxcv, maxcv_stop = maxcv_stop, 
            maxcv_autocheck = TRUE, jam.nM.iter = jam.nM.iter,extra.java.arguments=extra.java.arguments)
        binout <- as.matrix(topmods[, -ncol(topmods)])
        colnames(binout) <- colnames(topmods)[-ncol(topmods)]
        #snpmods <- apply(binout, 1,flashfm:::mod.fn)
        snpmods <- apply(binout, 1,mod.fn)
        nmod <- apply(binout, 1, sum)
        PP <- topmods[, ncol(topmods)]
        snpPP <- data.frame(rank = 1:length(nmod), size = nmod, 
            logPP = log(PP), PP = PP, str = snpmods, snps = snpmods, 
            stringsAsFactors = FALSE)
        snpPP <- snpPP[order(snpPP$PP, decreasing = TRUE), ]
        #expmods <- rlist::list.stack(lapply(snpPP$snps, flashfm:::tagexpand.mod,taglist = taglist))
        expmods <- rlist::list.stack(lapply(snpPP$snps, tagexpand.mod,taglist = taglist))
        wh <- which(duplicated(expmods$snps))
        if (length(wh) > 0) {
            expmods <- expmods[-wh, ]
        }
        row.names(expmods) <- expmods$snps
        check <- sapply(strsplit(expmods[, 2], "%"), function(x) length(x) > 
            length(unique(x)))
        if (sum(check) > 0) 
            expmods <- expmods[-which(check), ]
        #mbeta <- lapply(expmods[, 2], flashfm:::multibeta, beta1[[j]], 
        mbeta <- lapply(expmods[, 2], multibeta, beta1[[j]],                  
            covX, N = N, ybar = ybar, is.snpmat = FALSE, 
            raf = raf)
        names(mbeta) <- expmods[, 2]
        SSy <- Vy * (N - 1) + N * ybar^2
        Vx <- diag(covX)
        Mx <- 2 * raf
        #Sxy <- c(flashfm:::Sxy.hat(beta1 = beta1[[j]], Mx = Mx, N = N, 
        Sxy <- c(Sxy.hat(beta1 = beta1[[j]], Mx = Mx, N = N,                             
            Vx = Vx, muY = ybar), `1` = ybar * N)
        names(Sxy)[length(Sxy)] <- "one"
        #lABF <- sapply(expmods$snps, flashfm:::calcABF, mbeta, SSy = SSy, 
        lABF <- sapply(expmods$snps, calcABF, mbeta, SSy = SSy, 
            Sxy = Sxy, Vy = Vy, N = N)
        names(lABF) <- expmods$snps
        wh <- which(expmods$snps == "1")
        if (is.null(wh)) {
            dd <- data.frame(model = c("1", expmods$snps), 
                tag = c(FALSE, expmods$tag), lBF = c(0, lABF), 
                stringsAsFactors = FALSE)
            #l1 <- flashfm:::multibeta("1", beta1[[j]], covX, N = N, 
            l1 <- multibeta("1", beta1[[j]], covX, N = N, 
                ybar = ybar, is.snpmat = FALSE, raf = raf)
            mbeta <- rlist::list.append(mbeta, `1` = l1)
        }
        else {
            dd <- data.frame(model = expmods$snps, tag = expmods$tag, 
                lBF = lABF, stringsAsFactors = FALSE)
        }
        #SM <- flashfm:::makesnpmod(dd, expected = 2, nsnps = nsnps)
        SM <- makesnpmod(dd, expected = 2, nsnps = nsnps)
        
        #SM <- flashfm:::best.models.cpp(SM, maxmod = 1000)[[1]]
        SM <- best.models.cpp(SM, maxmod = 1000)[[1]]
        
        #SM <- flashfm:::PP2snpmod(SM)
        SM <- PP2snpmod(SM)
 
 		return(list(mbeta=mbeta, SM=SM))
}


#calc.maxmin <- MGflashfm:::calc.maxmin



#' @title Wrapper for JAMdynamic single-trait fine-mapping that also outputs SNP groups; PP and credible sets include SNP group information
#' @param gwas data.frame with 3 columns named: rsID, beta, EAF; 
#' @param N scalar that is the trait sample size 
#' @param corX SNP correlation matrix that must have snp names in row.names and col.names
#' @param save.path Path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param cred probability for credible set; default is 0.99
#' @param jam.nM.iter in millions, number of iterations to use in JAM; default 1 (1 million)
#' @param maxcv starting value for maximum number of causal variants
#' @param maxcv_stop maximum value to consider for maximum number of causal variants; maxcv_stop >= maxcv.
#' @param min.mppi trim snp groups with total MPPI < min.mppi in all diseases; default 0.01
#' @param minsnpmppi only group snps with total MPPI > minsnpmppi; default 0.01
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.5
#' @param extra.java.arguments A character string to be passed through to the java command line. E.g. to specify a
#' different temporary directory by passing "-Djava.io.tmpdir=/Temp".
#' @return List consisting of three objects: CS, a data.frame detailing the SNPs in the credible set; mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups, 
#' a list giving the SNP groups construced under single-trait fine-mapping
#' @export
#' @author Jenn Asimit
JAMdwithGroups <- function(gwas,N,corX, save.path, cred = 0.99, jam.nM.iter =1, maxcv = 1, maxcv_stop = 20,min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.6,extra.java.arguments =NULL) {

 fm <- JAMdynamic(gwas=gwas, corX=corX, ybar = 0, Vy = 1, N=N, cred=cred, save.path=save.path, 
    maxcv, maxcv_stop, jam.nM.iter,extra.java.arguments =extra.java.arguments) 
 ppdf <- fm$PP
 ppdf$str <- rownames(ppdf)   
 snpgroups <- makeSNPgroupsST(ppdf,corX, min.mppi, minsnpmppi, r2.minmerge) 
 mpp.pp <- PPsummariseST(ppdf,snpgroups)
 
 return(list(CS=data.frame(MPP=fm$CS,snp=names(fm$CS),
                           #group=flashfm::groupIDs.fn(snpgroups,names(fm$CS))),
                           group=groupIDs.fn(snpgroups,names(fm$CS))),
             mpp.pp=mpp.pp,snpGroups=snpgroups))
}



makeSNPgroupsST <- function (ppdf,corX, min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.5) 
{
    
    Xmat <- as.matrix(corX)
    sg <- groupU(ppdf, Xmat, is.snpmat=FALSE, min.mppi, minsnpmppi = minsnpmppi, 
        r2.minmerge)
    if (length(sg) > 1) {
        snpgroups <- sg$groups@.Data
        ng <- length(snpgroups)
        names(snpgroups) <- LETTERS[1:ng]
        if (ng > 26) 
            names(snpgroups)[27:min(ng, 52)] <- paste0(LETTERS[1:min(26, 
                ng - 26)], 2)
        if (ng > 52) {
            names(snpgroups)[27:min(ng, 52)] <- paste0(LETTERS[1:min(26, 
                ng - 26)], 2)
            names(snpgroups)[53:min(ng, 78)] <- paste0(LETTERS[1:min(26, 
                ng - 52)], 3)
        }
        if (ng > 78) {
            names(snpgroups)[27:min(ng, 52)] <- paste0(LETTERS[1:min(26, 
                ng - 26)], 2)
            names(snpgroups)[53:min(ng, 78)] <- paste0(LETTERS[1:min(26, 
                ng - 52)], 3)
            names(snpgroups)[79:min(ng, 104)] <- paste0(LETTERS[1:min(26, 
                ng - 78)], 4)
        }
        Ng <- lapply(snpgroups, length)
        sgd <- t(data.frame(Ng))
        colnames(sgd) <- "Group Size"
        message("SNP group sizes are: ")
        print(sgd)
    }
    else {
        snpgroups <- sg
    }
    return(snpgroups)
}



groupU <- function (ppdf, snp.data, is.snpmat, min.mppi = 0.01, minsnpmppi = 0.01, 
    r2.minmerge = 0.6) 
{
    #SM2 <- flashfm:::PP2snpmod(ppdf)
    SM2 <- PP2snpmod(ppdf)
    nsnps <- ncol(snp.data)
    s <- SM2@snps$var
    snps = colnames(snp.data)
    es <- setdiff(snps, s)
    if (length(es) > 0) {
            esmods <- data.frame(str = es, PP = 0, stringsAsFactors = FALSE)
            dd <- rbind(SM2@models[, c("str", "PP")], esmods)
            rownames(dd) <- dd$str
            #SM2 <- flashfm:::PP2snpmod(dd)
            SM2 <- PP2snpmod(dd)
            }
            
    #bs <- flashfm:::best.snps(SM2, pp.thr = 0)
    bs <- best.snps(SM2, pp.thr = 0)
    snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), "1")
    if (length(snps) == 1) {
        snpGroups <- list(A = snps)
        out <- snpGroups
    }
    if (length(snps) == 0) {
        minsnpmppi = 0
        min.mppi = 0
        snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, 
            ]$var), "1")
    }
    if (length(snps) == 1) {
        snpGroups <- list(A = snps)
        out <- snpGroups
    }
    if (length(snps) > 1) {
        if (is.snpmat) {
            snp.data <- snp.data[, snps]
            r2 <- cor(snp.data)^2
        }
        else {
            snp.data <- snp.data[snps, snps]
            r2 <- cov2cor(as.matrix(snp.data))^2
        }
        s <- rownames(dd)
        es <- setdiff(snps, s)
        #X <- flashfm:::makex(SM2)
        X <- makex(SM2)
        #mppi <- flashfm:::makemppi(X)
        mppi <- makemppi(X)
        MPPI <- matrix(mppi, nrow = 1)
        colnames(MPPI) <- names(mppi)
        
        MPPI <- t(MPPI)
        #R <- flashfm:::maker(X)[snps, snps]
        R <- maker(X)[snps, snps]
        rmax <- rmin <- R
#        if (length(R) > 1) 
#            for (i in 2:length(R)) rmin <- pmin(rmin, R[[i]])
#        rmax <- R[[1]]
#        if (length(R) > 1) 
#            for (i in 2:length(R)) rmax <- pmax(rmax, R[[i]])
        r <- ifelse(abs(rmin) > abs(rmax), rmin, rmax)
        #rd <- flashfm:::maked(R, r2)
        rd <- maked(R, r2)
        h <- hclust(rd, method = "complete")
        d <- as.dendrogram(h)
        r.tol = max(c(quantile(R, 0.9), 0))
        mem.sum <- function(members) {
            (sum(MPPI[members, , drop = FALSE]))
        }
        mem.marg <- function(members) {
                sum(pmin(apply(X$x[, members, drop = FALSE], 
                  1, sum), 1) * X$w)
            
        }
        mem.ab <- function(members) {
            list(a = mem.sum(members), b = mem.marg(members))
        }
        obj.ab <- function(object) {
            members <- labels(object)
            mem.ab(members)
        }
        mem.maxr.minr2 <- function(members) {
            r.sub <- r[members, members, drop = FALSE]
            r2.sub <- r2[members, members, drop = FALSE]
            mx <- max(r.sub[lower.tri(r.sub)], na.rm = TRUE)
            mn <- min(r2.sub[upper.tri(r2.sub)], na.rm = TRUE)
            c(mx, mn)
        }
        cutter <- function(object, mppi.max = 1.01, max.size = 50, 
            marg.sum.ratio = 1.1, max.r = 0, min.r2 = 0.5) {
            if (is.leaf(object)) 
                return(labels(object))
            members <- labels(object)
            ab <- mem.ab(members)
            if (max(ab[[1]]) < min.mppi) 
                return(labels(object))
            if (max(ab[[1]]) > mppi.max) 
                return(list(cutter(object[[1]], min.r2 = min.r2), 
                  cutter(object[[2]], min.r2 = min.r2)))
            mxmn <- mem.maxr.minr2(members)
            if (mxmn[1] > r.tol || mxmn[2] < min.r2) 
                return(list(cutter(object[[1]], min.r2 = min.r2), 
                  cutter(object[[2]], min.r2 = min.r2)))
            if (min(c(mem.sum(labels(object[[1]])), mem.sum(labels(object[[2]]))) < 
                min.mppi)) 
                return(list(cutter(object[[1]], min.r2 = min.r2), 
                  cutter(object[[2]], min.r2 = min.r2)))
            if (max(ab[[1]]) <= mppi.max & all(ab[[1]] < ab[[2]] * 
                marg.sum.ratio)) 
                return(labels(object))
            return(list(cutter(object[[1]], min.r2 = min.r2), 
                cutter(object[[2]], min.r2 = min.r2)))
        }
        mem.summ <- function(members) {
            n <- length(members)
            ab <- mem.ab(members)
            mppi.min <- apply(MPPI[members, , drop = FALSE], 
                2, min)
            mppi.max <- apply(MPPI[members, , drop = FALSE], 
                2, max)
            r2.sub <- r2[members, members]
            r2.summ <- summary(r2.sub[upper.tri(r2.sub)])
            r.sub <- r[members, members]
            r.summ <- summary(r.sub[upper.tri(r.sub)])
            c(n = n, sum.mppi = ab[[1]], r2 = r2.summ["Min."], 
                r2 = r2.summ["Max."], r = r.summ["Min."], r = r.summ["Max."], 
                mppi.min = mppi.min, mppi.max = mppi.max)
        }
        ret <- cutter(d, min.r2 = r2.minmerge)
        if (!is.list(ret)) 
            ret <- list(ret)
        #ret <- flashfm:::LinearizeNestedList(ret)
        ret <- LinearizeNestedList(ret)
        ret.mppi <- t(sapply(ret, mem.sum))
        use <- apply(ret.mppi, 1, max) > minsnpmppi
        df <- sapply(ret, mem.summ)
        df <- t(df)
        union.summary <- df[use, , drop = FALSE]
        union.content <- ret[use]
        use <- apply(union.summary[, grep("sum.mppi", colnames(union.summary)), 
            drop = FALSE], 1, max) > min.mppi
        G1 <- union.summary[use, , drop = FALSE]
        G2 <- union.content[use]
        rownames(G1) <- NULL
        merger <- function(G1, G2) {
            maxr2 <- calc.maxmin(r2, G2, fun = max)
            minr2 <- calc.maxmin(r2, G2, fun = min)
            maxr <- calc.maxmin(r, G2, fun = max)
            diag(maxr2) <- 0
            tomerge <- maxr2 > r2.minmerge & maxr < r.tol
            if (any(tomerge, na.rm = TRUE)) {
                wh <- which(tomerge, arr.ind = TRUE)
                wh <- wh[wh[, 1] < wh[, 2], , drop = FALSE]
                wh <- cbind(wh, maxr2[wh])
                wh <- wh[order(wh[, 3], decreasing = TRUE), , 
                  drop = FALSE]
                for (k in 1:nrow(wh)) {
                  a <- wh[k, 1]
                  b <- wh[k, 2]
                  sumcols <- grep("sum.mppi", colnames(G1))
                  if (any(colSums(G1[c(a, b), sumcols, drop = FALSE]) > 
                    1.01)) 
                    next
                  G2[[a]] <- c(G2[[a]], G2[[b]])
                  G2[[b]] <- NULL
                  for (nm in c(1, sumcols)) G1[a, nm] <- sum(G1[c(a, 
                    b), nm])
                  for (nm in setdiff(1:ncol(G1), c(1, sumcols))) G1[a, 
                    nm] <- max(G1[c(a, b), nm])
                  G1 <- G1[-b, , drop = FALSE]
                  return(merger(G1, G2))
                }
            }
            return(list(G1, G2))
        }
        G <- merger(G1, G2)
        tmp <- G[[2]]
        names(tmp) <- sapply(tmp, "[[", 1)
        newgroups <- new("groups", tmp, tags = names(tmp))
        out <- list(summary = G[[1]], groups = newgroups, r2 = r2)
    }
    return(out)
}





PPsummariseST <- function (ppdf, snpgroups) 
{
    pp <- ppdf[,"PP",drop=F]
    #mpp <- MGflashfm:::MPPcalc(pp)
    mpp <- MPPcalc(pp)
  
    ppvec <- pp[,1]
    names(ppvec) <- rownames(pp)
    #PPout <- flashfm:::PPmodGroups(ppvec, snpgroups, minPP = 0)
    PPout <- PPmodGroups(ppvec, snpgroups, minPP = 0)
    PPg <- PPout$group[,"PP",drop=F]
    #MPPg <- MGflashfm:::MPPcalc(PPg)
    MPPg <- MPPcalc(PPg)
    
    mpp <- data.frame(mpp[order(mpp,decreasing=T)])
    MPPg <- data.frame(MPPg[order(MPPg,decreasing=T)])
    colnames(mpp)=colnames(MPPg) = "MPP"
       
    return(list(MPP = mpp, MPPg = MPPg, PP = pp, PPg = PPg))
}


