
#' @title Make SNP groups using fine-mapping information from all of the traits
#' @param main.input output from flashfm.input function
#' @param is.snpmat logical taking value TRUE when genotype matrix is provided and FALSE when covariance matrix is given
#' @param min.mppi trim snp groups with total MPPI < min.mppi in all diseases; default 0.01
#' @param minsnpmppi min snp mppi
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.5
#' @return list of  SNP groups
#' @export
#copied from makeSNPgroupsU
makeSNPgroups <- function (main.input, is.snpmat, min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.5) 
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
#previous version
# makeSNPgroups <- function(main.input,is.snpmat,min.mppi = 0.01,r2.minmerge=0.5) {
# snp.data <- main.input$Gmat
# SMlist <- main.input$SM
# #if(is.snpmat) { Xmat <- new("SnpMatrix",round(snp.data+1)) 
# #} else {Xmat <- as.matrix(snp.data) } 
# Xmat <- as.matrix(snp.data)
# sg <- groupmulti(SMlist,Xmat,is.snpmat,min.mppi,minsnpmppi=.001,r2.minmerge)
# 
# if(length(sg) > 1){
# snpgroups <- sg$groups@.Data
# ng <- length(snpgroups)
# names(snpgroups) <- LETTERS[1:ng] # arbitrary names
# if(ng>26) names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
# if(ng>52) { 
#  names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
#  names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
#  }
# if(ng>78) {
#  names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
#  names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
#  names(snpgroups)[79:min(ng,104)] <- paste0(LETTERS[1:min(26,ng-78)],4)
# }
# Ng <- lapply(snpgroups,length)
# sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
# message("SNP group sizes are: "); print(sgd)
# } else { snpgroups <- sg}
# return(snpgroups)
# }


#' @title Make two sets of SNP groups using fine-mapping information from all of the traits using two sets of results and maps the names between them 
#' @param main.input output from flashfm.input function
#' @param fm.multi output from flashfm function
#' @param is.snpmat logical taking value TRUE when genotype matrix is provided and FALSE when covariance matrix is given
#' @param min.mppi trim snp groups with total MPPI < min.mppi in all diseases; default 0.01
#' @param minsnpmppi only group snps with total MPPI > minsnpmppi; default 0.001
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.5
#' @return list of  three objects: groups.fm is a list of SNP groups using the single-trait results; groups.flashfm is a list of SNP groups using the flashfm results; group.sizes is  a table of SNP group sizes for the two sets of groups
#' @export
#copied from makeSNPgroups2U
makeSNPgroups2 <- function(main.input,fm.multi,is.snpmat,min.mppi = 0.01,minsnpmppi=0.001,r2.minmerge=0.5) {
  snp.data <- main.input$Gmat
  M <- length(fm.multi$PP)
  #SMlist <- main.input$SM
  #if(is.snpmat) { Xmat <- new("SnpMatrix",round(snp.data+1)) 
  #} else {Xmat <- as.matrix(snp.data) } 
  SMlist <- vector("list",M)
  fmpp  <- fm.multi$PP
  for(i in 1:M) {
    ppdf <- data.frame(str=as.character(rownames(fmpp[[i]])),PP=fmpp[[i]][,1], stringsAsFactors = FALSE)
    SMlist[[i]] <- PP2snpmod(ppdf)
  }
  Xmat <- as.matrix(snp.data)
  sg <- groupmultiU(SMlist,Xmat,is.snpmat,min.mppi,minsnpmppi,r2.minmerge)
  if(length(sg) > 1){
    snpgroups <- sg$groups@.Data
    ng <- length(snpgroups)
    names(snpgroups) <- LETTERS[1:ng] # arbitrary names
    if(ng>26) names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:(ng-26)],2)
    if(ng>52) { 
      names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
      names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
    }
    if(ng>78) {
      names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
      names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
      names(snpgroups)[79:min(ng,104)] <- paste0(LETTERS[1:min(26,ng-78)],4)
    }
    Ng <- lapply(snpgroups,length)
    #sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
    #message("SNP group sizes based on single-trait results are: "); print(sgd)
  } else { 
    snpgroups <- sg
    ng <- 1
  }
  
  fmpp  <- fm.multi$PP
  M <- length(fmpp)
  SM2list <- vector("list",M)
  for(i in 1:M) {
    ppdf <- data.frame(str=as.character(rownames(fmpp[[i]])),PP=fmpp[[i]][,2], stringsAsFactors = FALSE)
    SM2list[[i]] <- PP2snpmod(ppdf)
  }
  sg2 <- groupmultiU(SM2list,Xmat,is.snpmat,min.mppi,minsnpmppi,r2.minmerge)
  if(length(sg2) > 1){
    snpgroups2 <- sg2$groups@.Data
    ng2 <- length(snpgroups2)
    names(snpgroups2) <- LETTERS[1:ng2] # arbitrary names
    if(ng2>26) names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
    if(ng2>52) { 
      names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
      names(snpgroups2)[53:min(ng2,78)] <- paste0(LETTERS[1:min(26,ng2-52)],3)
    }
    if(ng2>78) {
      names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
      names(snpgroups2)[53:min(ng2,78)] <- paste0(LETTERS[1:min(26,ng2-52)],3)
      names(snpgroups2)[79:min(ng2,104)] <- paste0(LETTERS[1:min(26,ng2-78)],4)
    }
    
    Ng2 <- lapply(snpgroups2,length)
    #sgd2 <- t(data.frame(Ng2)); colnames(sgd2) <- "Group Size"
    #message("SNP group sizes based on flashfm results are: "); print(sgd2)
  } else { 
    snpgroups2 <- sg2
    ng2 <- 1
  }
  
  wh <- NULL
  for(i in 1:ng) {
    for(j in 1:ng2) {
      if(length(intersect(snpgroups[[i]],snpgroups2[[j]])) > 0 )wh <- rbind(wh,c(i,j))
    }
  }
  
  
  
  newsg2 <- snpgroups2  
  newnames <- character(ng2)
  
  ind <- which(1:ng %in% wh[,1])
  wrm <- c() 
  
  if(any(duplicated(wh[,1]))) {
    dup <- which(duplicated(wh[,1]))
    dups <- c()
    
    for(i in dup) { 
      dd <- which(wh[,1]==wh[i,1])
      for(j in 1:length(dd)) newnames[wh[dd[j],2]] <- paste(names(snpgroups)[wh[dd[j],1]],j,sep=".")  			
      wrm <- c(wrm,dd)
    }
    wh <- matrix(wh[-wrm,],ncol=2)
  }
  
  
  for(i in 1:nrow(wh)) newnames[wh[i,2]] <- names(snpgroups)[wh[i,1]]
  
  
  
  nc <- nchar(newnames)
  snames <- c(LETTERS[1:26],paste0(LETTERS[1:26],2),paste0(LETTERS[1:26],3),paste0(LETTERS[1:26],4))
  if(any(nc==0)){
    ind <- which(nc==0)
    newnames[ind] <-  snames[(ng+1):(ng+length(ind))]	
  }     
  names(snpgroups2) <- newnames 
  snpgroups2 <- snpgroups2[order(names(snpgroups2))]
  
  Ng <- lapply(snpgroups,length)
  sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
  message("SNP group sizes based on single-trait results are: "); print(sgd)
  
  Ng2 <- lapply(snpgroups2,length)
  sgd2 <- t(data.frame(Ng2)); colnames(sgd2) <- "Group Size"
  message("SNP group sizes based on flashfm results are: "); print(sgd2)
  
  group.sizes <- do.call("smartbind",c(list(t(sgd),t(sgd2)),fill=0))
  rownames(group.sizes) <- fm.multi$sharing
  
  return(list(groups.fm=snpgroups,groups.flashfm=snpgroups2, group.sizes=group.sizes))
}
#previous version
# makeSNPgroups2 <- function(main.input,fm.multi,is.snpmat,min.mppi = 0.01,minsnpmppi=0.001,r2.minmerge=0.5) {
# snp.data <- main.input$Gmat
# 
# #SMlist <- main.input$SM
# #if(is.snpmat) { Xmat <- new("SnpMatrix",round(snp.data+1)) 
# #} else {Xmat <- as.matrix(snp.data) } 
# M <- length(fm.multi$PP)
# SMlist <- vector("list",M)
# fmpp  <- fm.multi$PP
# for(i in 1:M) {
#  ppdf <- data.frame(str=as.character(rownames(fmpp[[i]])),PP=fmpp[[i]][,1], stringsAsFactors = FALSE)
#  SMlist[[i]] <- PP2snpmod(ppdf)
#  }
# Xmat <- as.matrix(snp.data)
# sg <- groupmulti(SMlist,Xmat,is.snpmat,min.mppi,minsnpmppi,r2.minmerge)
# if(length(sg) > 1){
# snpgroups <- sg$groups@.Data
# ng <- length(snpgroups)
# names(snpgroups) <- LETTERS[1:ng] # arbitrary names
# if(ng>26) names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:(ng-26)],2)
# if(ng>52) { 
#  names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
#  names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
#  }
# if(ng>78) {
#  names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
#  names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
#  names(snpgroups)[79:min(ng,104)] <- paste0(LETTERS[1:min(26,ng-78)],4)
# }
# Ng <- lapply(snpgroups,length)
# #sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
# #message("SNP group sizes based on single-trait results are: "); print(sgd)
# } else { 
#  snpgroups <- sg
#  ng <- 1
#  }
# 
# fmpp  <- fm.multi$PP
# M <- length(fmpp)
# SM2list <- vector("list",M)
# for(i in 1:M) {
#  ppdf <- data.frame(str=as.character(rownames(fmpp[[i]])),PP=fmpp[[i]][,2], stringsAsFactors = FALSE)
#  SM2list[[i]] <- PP2snpmod(ppdf)
#  }
# sg2 <- groupmulti(SM2list,Xmat,is.snpmat,min.mppi,minsnpmppi,r2.minmerge)
# if(length(sg2) > 1){
# snpgroups2 <- sg2$groups@.Data
# ng2 <- length(snpgroups2)
# names(snpgroups2) <- LETTERS[1:ng2] # arbitrary names
# if(ng2>26) names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
# if(ng2>52) { 
#  names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
#  names(snpgroups2)[53:min(ng2,78)] <- paste0(LETTERS[1:min(26,ng2-52)],3)
#  }
# if(ng2>78) {
#  names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
#  names(snpgroups2)[53:min(ng2,78)] <- paste0(LETTERS[1:min(26,ng2-52)],3)
#  names(snpgroups2)[79:min(ng2,104)] <- paste0(LETTERS[1:min(26,ng2-78)],4)
# }
# 
# Ng2 <- lapply(snpgroups2,length)
# #sgd2 <- t(data.frame(Ng2)); colnames(sgd2) <- "Group Size"
# #message("SNP group sizes based on flashfm results are: "); print(sgd2)
# } else { 
#  snpgroups2 <- sg2
#  ng2 <- 1
#  }
# 
# wh <- NULL
# for(i in 1:ng) {
#  for(j in 1:ng2) {
#   if(length(intersect(snpgroups[[i]],snpgroups2[[j]])) > 0 )wh <- rbind(wh,c(i,j))
#  }
# }
# 
# 
#  
#  newsg2 <- snpgroups2  
#  newnames <- character(ng2)
#  
#  ind <- which(1:ng %in% wh[,1])
#  wrm <- c() 
# 
#  if(any(duplicated(wh[,1]))) {
#   dup <- which(duplicated(wh[,1]))
#   dups <- c()
#   
#   for(i in dup) { 
#   	dd <- which(wh[,1]==wh[i,1])
#   	for(j in 1:length(dd)) newnames[wh[dd[j],2]] <- paste(names(snpgroups)[wh[dd[j],1]],j,sep=".")  			
# 	wrm <- c(wrm,dd)
# 	}
#    wh <- matrix(wh[-wrm,],ncol=2)
#    }
#    
#    
# for(i in 1:nrow(wh)) newnames[wh[i,2]] <- names(snpgroups)[wh[i,1]]
#  
# 
# 
# nc <- nchar(newnames)
# snames <- c(LETTERS[1:26],paste0(LETTERS[1:26],2),paste0(LETTERS[1:26],3),paste0(LETTERS[1:26],4))
# if(any(nc==0)){
# 	ind <- which(nc==0)
# 	newnames[ind] <-  snames[(ng+1):(ng+length(ind))]	
# 	}     
# names(snpgroups2) <- newnames 
# snpgroups2 <- snpgroups2[order(names(snpgroups2))]
# 
# Ng <- lapply(snpgroups,length)
# sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
# message("SNP group sizes based on single-trait results are: "); print(sgd)
# 
# Ng2 <- lapply(snpgroups2,length)
# sgd2 <- t(data.frame(Ng2)); colnames(sgd2) <- "Group Size"
# message("SNP group sizes based on flashfm results are: "); print(sgd2)
# 
# group.sizes <- do.call("smartbind",c(list(t(sgd),t(sgd2)),fill=0))
# rownames(group.sizes) <- fm.multi$sharing
# 
# return(list(groups.fm=snpgroups,groups.flashfm=snpgroups2, group.sizes=group.sizes))
# }



PP2snpmod <- function (ppdf) 
# modified abf2snpmod from Chris Wallace
# have data.frame with columns "str" = models with SNPs separated by %; "PP" model PP
{
    tmp <- new("snpmod")
    msize <- nchar(gsub("[^%]", "", ppdf$str)) + 1
    msize[ppdf$str == "1"] <- 0
    
       tmp@models <- data.frame(str = ppdf$str, size = msize, 
        lPP = log(ppdf$PP), PP = ppdf$PP, stringsAsFactors = FALSE)
    tmp@model.snps <- strsplit(tmp@models$str, "%")
    message("calculating marginal SNP inclusion probabilities")
    marg.snps(tmp)
}



##' @title Group SNPs; adapted from group.multi by Chris Wallace
##' @param SM2 snpmod object
##' @param snp.data SnpMatrix or SNP covariance matrix with snp ids as column names 
##' @param is.snpmat logical taking value TRUE when SnpMatrix is provided and FALSE when covariance matrix is given
##' @param min.mppi trim snp groups with total MPPI < min.mppi in all
##'     diseases
##' @param minsnpmppi only group snps with total MPPI > minsnpmppi
##' @param r2.minmerge merge groups with minimum between-group r2 >
##'     r2.minmerge
##' @return list with three components.
##'
##' First is a data.frame with each row giving summary statistics for
##'     each group.
##'
##' Second is a groups object, each elements ordered according to the rows of the summary
##' 
##' Third is the r2 matrix calculated.
#' @export
#copied from groupmultiU
groupmulti <- function (SM2, snp.data, is.snpmat, min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.6) 
{
  stopifnot(is.list(SM2))
  M <- length(SM2)
  dd <- vector("list",M)
  
  if(is(SM2[[1]], "snpmod")) {
    nsnps <- ncol(snp.data)
    s <- lapply(SM2,function(x) x@snps$var)
    es <- lapply(s,function(s,snps) setdiff(snps,s),snps=colnames(snp.data))
    for(j in 1:M) {
      if(length(es[[j]])>0){
        esmods <- data.frame(str=es[[j]],PP=0,stringsAsFactors = FALSE)
        dd[[j]] <- rbind(SM2[[j]]@models[,c("str","PP")],esmods)
        rownames(dd[[j]]) <- dd[[j]]$str
        SM2[[j]] <- PP2snpmod(dd[[j]])
      }
      
    }
    
  }
  
  
  bs <- best.snps(SM2, pp.thr = 0)
  bs <- do.call("rbind", bs)
  
  snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), 
                  "1")
  #    if(length(snps) < 2) {
  #     minsnpmppi = 0
  #     min.mppi = 0
  #     snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), 
  #        "1")
  #     }
  if(length(snps) == 1) {
    snpGroups <- list(A=snps)
    out <- snpGroups
  }
  if(length(snps) == 0) {
    minsnpmppi = 0
    min.mppi = 0
    snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), "1")
  }
  if(length(snps) == 1) {
    snpGroups <- list(A=snps)
    out <- snpGroups
  }
  
  
  if(length(snps) > 1) {    
    if(is.snpmat) {
      snp.data <- snp.data[, snps]
      r2 <- cor(snp.data)^2
    } else { 
      snp.data <- snp.data[snps,snps]
      r2 <- cov2cor(as.matrix(snp.data))^2 
    }
    
    #	s <- lapply(SM2,function(x) x@snps$var)
    s <- lapply(dd, function(x) rownames(dd))
    es <- lapply(s,function(s,snps) setdiff(snps,s),snps=snps) # prioritised snps not among models for a trait
    
    
    X <- lapply(SM2, makex)
    mppi <- lapply(X, makemppi)
    
    MPPI <- do.call("smartbind",c(lapply(mppi,function(m) {mm <- matrix(m,nrow=1); colnames(mm) <- names(m); return(mm)}),fill=0))
    MPPI <- t(MPPI)
    ##    MPPI <- do.call("cbind", lapply(X, makemppi))
    R <- lapply(X, function(x) maker(x)[snps, snps])
    rmax <- rmin <- R[[1]]
    if (length(R) > 1) 
      for (i in 2:length(R)) rmin <- pmin(rmin, R[[i]])
    rmax <- R[[1]]
    if (length(R) > 1) 
      for (i in 2:length(R)) rmax <- pmax(rmax, R[[i]])
    r <- ifelse(abs(rmin) > abs(rmax), rmin, rmax)
    rd <- maked(r, r2)
    h <- hclust(rd, method = "complete")
    d <- as.dendrogram(h)
    r.tol = max(c(quantile(r, 0.9),0))
    
    ## utility functions in local environment
    mem.sum <- function(members) {
      (colSums(MPPI[members, , drop = FALSE]))
    }
    mem.marg <- function(members) {
      sapply(X, function(x) {
        sum(pmin(apply(x$x[, members, drop = FALSE], 1, sum), 
                 1) * x$w)
      })
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
        return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
      mxmn <- mem.maxr.minr2(members)
      if (mxmn[1] > r.tol || mxmn[2] < min.r2) 
        return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
      if (min(c(mem.sum(labels(object[[1]])), mem.sum(labels(object[[2]]))) < 
              min.mppi)) 
        return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
      if (max(ab[[1]]) <= mppi.max & all(ab[[1]] < ab[[2]] * 
                                         marg.sum.ratio)) 
        return(labels(object))
      return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
    }
    mem.summ <- function(members) {
      n <- length(members)
      ab <- mem.ab(members)
      mppi.min <- apply(MPPI[members, , drop = FALSE], 2, min)
      mppi.max <- apply(MPPI[members, , drop = FALSE], 2, max)
      r2.sub <- r2[members, members]
      r2.summ <- summary(r2.sub[upper.tri(r2.sub)])
      r.sub <- r[members, members]
      r.summ <- summary(r.sub[upper.tri(r.sub)])
      c(n = n, sum.mppi = ab[[1]], r2 = r2.summ["Min."], r2 = r2.summ["Max."], 
        r = r.summ["Min."], r = r.summ["Max."], mppi.min = mppi.min, 
        mppi.max = mppi.max)
    }
    ret <- cutter(d,min.r2=r2.minmerge)
    if (!is.list(ret)) 
      ret <- list(ret)
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
        wh <- wh[order(wh[, 3], decreasing = TRUE), , drop = FALSE]
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
#previous version
# groupmulti <- function (SM2, snp.data, is.snpmat, min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.6) 
# {
#     stopifnot(is.list(SM2))
#     M <- length(SM2)
#     dd <- vector("list",M)
# 	
#     if(is(SM2[[1]], "snpmod")) {
#       nsnps <- ncol(snp.data)
#       s <- lapply(SM2,function(x) x@snps$var)
# 	es <- lapply(s,function(s,snps) setdiff(snps,s),snps=colnames(snp.data))
# 	for(j in 1:M) {
# 		if(length(es[[j]])>0){
# 		esmods <- data.frame(str=es[[j]],PP=0,stringsAsFactors = FALSE)
# 		dd[[j]] <- rbind(SM2[[j]]@models[,c("str","PP")],esmods)
# 		rownames(dd[[j]]) <- dd[[j]]$str
# 		SM2[[j]] <- PP2snpmod(dd[[j]])
#  				}
#  				
# 		}
#     
#     }
#        
#     
#     bs <- best.snps(SM2, pp.thr = 0)
#     bs <- do.call("rbind", bs)
# 
#     snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), 
#         "1")
# #    if(length(snps) < 2) {
# #     minsnpmppi = 0
# #     min.mppi = 0
# #     snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), 
# #        "1")
# #     }
#      if(length(snps) == 1) {
#       snpGroups <- list(A=snps)
#       out <- snpGroups
#      }
#      if(length(snps) == 0) {
#       minsnpmppi = 0
#       min.mppi = 0
#       snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), "1")
#      }
#      if(length(snps) == 1) {
#       snpGroups <- list(A=snps)
#       out <- snpGroups
#      }
#      
#      
#     if(length(snps) > 1) {    
#     if(is.snpmat) {
#     	snp.data <- snp.data[, snps]
# 		r2 <- cor(snp.data)^2
# 		} else { 
# 			snp.data <- snp.data[snps,snps]
# 			r2 <- cov2cor(as.matrix(snp.data))^2 
# 			}
# 
# #	s <- lapply(SM2,function(x) x@snps$var)
# 	s <- lapply(dd, function(x) rownames(dd))
# 	es <- lapply(s,function(s,snps) setdiff(snps,s),snps=snps) # prioritised snps not among models for a trait
# 	
# 
#     X <- lapply(SM2, makex)
#     mppi <- lapply(X, makemppi)
#     
#     MPPI <- do.call("smartbind",c(lapply(mppi,function(m) {mm <- matrix(m,nrow=1); colnames(mm) <- names(m); return(mm)}),fill=0))
#     MPPI <- t(MPPI)
# ##    MPPI <- do.call("cbind", lapply(X, makemppi))
#     R <- lapply(X, function(x) maker(x)[snps, snps])
#     rmax <- rmin <- R[[1]]
#     if (length(R) > 1) 
#         for (i in 2:length(R)) rmin <- pmin(rmin, R[[i]])
#     rmax <- R[[1]]
#     if (length(R) > 1) 
#         for (i in 2:length(R)) rmax <- pmax(rmax, R[[i]])
#     r <- ifelse(abs(rmin) > abs(rmax), rmin, rmax)
#     rd <- maked(r, r2)
#     h <- hclust(rd, method = "complete")
#     d <- as.dendrogram(h)
#     r.tol = max(c(quantile(r, 0.9),0))
#     
#     ## utility functions in local environment
#     mem.sum <- function(members) {
#         (colSums(MPPI[members, , drop = FALSE]))
#     }
#     mem.marg <- function(members) {
#         sapply(X, function(x) {
#             sum(pmin(apply(x$x[, members, drop = FALSE], 1, sum), 
#                 1) * x$w)
#         })
#     }
#     mem.ab <- function(members) {
#         list(a = mem.sum(members), b = mem.marg(members))
#     }
#     obj.ab <- function(object) {
#         members <- labels(object)
#         mem.ab(members)
#     }
#     mem.maxr.minr2 <- function(members) {
#         r.sub <- r[members, members, drop = FALSE]
#         r2.sub <- r2[members, members, drop = FALSE]
#         mx <- max(r.sub[lower.tri(r.sub)], na.rm = TRUE)
#         mn <- min(r2.sub[upper.tri(r2.sub)], na.rm = TRUE)
#         c(mx, mn)
#     }
#     cutter <- function(object, mppi.max = 1.01, max.size = 50, 
#         marg.sum.ratio = 1.1, max.r = 0, min.r2 = 0.5) {
#         if (is.leaf(object)) 
#             return(labels(object))
#         members <- labels(object)
#         ab <- mem.ab(members)
#         if (max(ab[[1]]) < min.mppi) 
#             return(labels(object))
#         if (max(ab[[1]]) > mppi.max) 
#             return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
#         mxmn <- mem.maxr.minr2(members)
#         if (mxmn[1] > r.tol || mxmn[2] < min.r2) 
#             return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
#         if (min(c(mem.sum(labels(object[[1]])), mem.sum(labels(object[[2]]))) < 
#             min.mppi)) 
#             return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
#         if (max(ab[[1]]) <= mppi.max & all(ab[[1]] < ab[[2]] * 
#             marg.sum.ratio)) 
#             return(labels(object))
#         return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
#     }
#     mem.summ <- function(members) {
#         n <- length(members)
#         ab <- mem.ab(members)
#         mppi.min <- apply(MPPI[members, , drop = FALSE], 2, min)
#         mppi.max <- apply(MPPI[members, , drop = FALSE], 2, max)
#         r2.sub <- r2[members, members]
#         r2.summ <- summary(r2.sub[upper.tri(r2.sub)])
#         r.sub <- r[members, members]
#         r.summ <- summary(r.sub[upper.tri(r.sub)])
#         c(n = n, sum.mppi = ab[[1]], r2 = r2.summ["Min."], r2 = r2.summ["Max."], 
#             r = r.summ["Min."], r = r.summ["Max."], mppi.min = mppi.min, 
#             mppi.max = mppi.max)
#     }
#     ret <- cutter(d,min.r2=r2.minmerge)
#     if (!is.list(ret)) 
#         ret <- list(ret)
#     ret <- LinearizeNestedList(ret)
#     ret.mppi <- t(sapply(ret, mem.sum))
#     use <- apply(ret.mppi, 1, max) > minsnpmppi
#     df <- sapply(ret, mem.summ)
#     df <- t(df)
#     union.summary <- df[use, , drop = FALSE]
#     union.content <- ret[use]
#     use <- apply(union.summary[, grep("sum.mppi", colnames(union.summary)), 
#         drop = FALSE], 1, max) > min.mppi
#     G1 <- union.summary[use, , drop = FALSE]
#     G2 <- union.content[use]
#     rownames(G1) <- NULL
#     merger <- function(G1, G2) {
#         maxr2 <- calc.maxmin(r2, G2, fun = max)
#         minr2 <- calc.maxmin(r2, G2, fun = min)
#         maxr <- calc.maxmin(r, G2, fun = max)
#         diag(maxr2) <- 0
#         tomerge <- maxr2 > r2.minmerge & maxr < r.tol
#         if (any(tomerge, na.rm = TRUE)) {
#             wh <- which(tomerge, arr.ind = TRUE)
#             wh <- wh[wh[, 1] < wh[, 2], , drop = FALSE]
#             wh <- cbind(wh, maxr2[wh])
#             wh <- wh[order(wh[, 3], decreasing = TRUE), , drop = FALSE]
#             for (k in 1:nrow(wh)) {
#                 a <- wh[k, 1]
#                 b <- wh[k, 2]
#                 sumcols <- grep("sum.mppi", colnames(G1))
#                 if (any(colSums(G1[c(a, b), sumcols, drop = FALSE]) > 
#                   1.01)) 
#                   next
#                 G2[[a]] <- c(G2[[a]], G2[[b]])
#                 G2[[b]] <- NULL
#                 for (nm in c(1, sumcols)) G1[a, nm] <- sum(G1[c(a, 
#                   b), nm])
#                 for (nm in setdiff(1:ncol(G1), c(1, sumcols))) G1[a, 
#                   nm] <- max(G1[c(a, b), nm])
#                 G1 <- G1[-b, , drop = FALSE]
#                 return(merger(G1, G2))
#             }
#         }
#         return(list(G1, G2))
#     }
#     G <- merger(G1, G2)
#     tmp <- G[[2]]
#     names(tmp) <- sapply(tmp, "[[", 1)
#     newgroups <- new("groups", tmp, tags = names(tmp))
#     out <- list(summary = G[[1]], groups = newgroups, r2 = r2)
#     }
#     
#     return(out)
# }

##

##

sparse.cor <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x)
    cSums <- colSums(x)
    ## Calculate the population covariance matrix.
    ## There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    ## The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(crossprod(x))
    covmat <- covmat+crossp
    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}

##' @importFrom Matrix sparseMatrix
makex <- function(obj) {
    snps <- sort(rownames(obj@snps))
    nums <- seq_along(snps)
    names(nums) <- snps
    npermod <- sapply(obj@model.snps,length)
    x <- sparseMatrix(i=rep(1:nrow(obj@models),times=npermod),
                      j=nums[unlist(obj@model.snps)],
                      x=rep(1,length(unlist(obj@model.snps))))
                                        #x=rep(obj@models$PP,times=npermod)^0.5)
    dimnames(x)[[2]] <- snps
    w <- obj@models$PP
    if(snps[[1]]=="1")
        x <- x[,-1]
    return(list(x=x,w=w))
}
makemppi <- function(x) {
    wx <- x$w*x$x
    colSums(wx)        
}
maker <- function(x) {
    sx <- sampx(x,n=1000000)
    r=sparse.cor(sx) #sqrt(w) * x) # sx
    diag(r) <-0 
    r[is.na(r)] <- 0
    r
}
maked <- function(r,r2) {
    ## rd <- dist(scale(cbind(r,r2)))
    ## rd <- dist(scale(r2))
    ## rd <- dist((sign(r)*(r2)))
    if(any(is.na(r2))) 
        r2[ which(is.na(r2)) ] <- 0
    if(any(r==0))
        r[ which(r==0) ] <- -0.001
    rd <- as.dist( ( sign(r) * (r2) + 1 ) / 2)
}
sampx <- function(d,n=1000) {
    wh <- sample(1:nrow(d$x),n,replace=TRUE,prob=d$w)
    d$x[wh,]
}

##' @title calculate max or min of subset of a matrix
##' @param r2 matrix, with row/column names which are members of L1, L2
##' @param L1 first list of vectors of snps
##' @param L2 second list vectors of snps
##' @param fun max or min, but could be another function that returns a scalar. unquoted.
##' @return smaller matrix indexed by elements of L1 (rows) and L2 (columns)
##' @export
##' @author Chris Wallace
calc.maxmin <- function(r2,L1,L2=L1,fun=max) {
    maxr2 <- matrix(0,length(L1),length(L2))
    for(i in seq_along(L1)) {
        for(j in seq_along(L2)) {
            maxr2[i,j] <- fun(r2[ L1[[i]], L2[[j]] ])
        }
    }
    maxr2
}

## https://gist.github.com/mrdwab/4205477
## LinearizeNestedList:
##
## https://sites.google.com/site/akhilsbehl/geekspace/
##         articles/r/linearize_nested_lists_in_r
##
## Akhil S Bhel
## 
## Implements a recursive algorithm to linearize nested lists upto any
## arbitrary level of nesting (limited by R's allowance for recursion-depth).
## By linearization, it is meant to bring all list branches emanating from
## any nth-nested trunk upto the top-level trunk s.t. the return value is a
## simple non-nested list having all branches emanating from this top-level
## branch.
##
## Since dataframes are essentially lists a boolean option is provided to
## switch on/off the linearization of dataframes. This has been found
## desirable in the author's experience.
##
## Also, one'd typically want to preserve names in the lists in a way as to
## clearly denote the association of any list element to it's nth-level
## history. As such we provide a clean and simple method of preserving names
## information of list elements. The names at any level of nesting are
## appended to the names of all preceding trunks using the `NameSep` option
## string as the seperator. The default `/` has been chosen to mimic the unix
## tradition of filesystem hierarchies. The default behavior works with
## existing names at any n-th level trunk, if found; otherwise, coerces simple
## numeric names corresponding to the position of a list element on the
## nth-trunk. Note, however, that this naming pattern does not ensure unique
## names for all elements in the resulting list. If the nested lists had
## non-unique names in a trunk the same would be reflected in the final list.
## Also, note that the function does not at all handle cases where `some`
## names are missing and some are not.
##
## Clearly, preserving the n-level hierarchy of branches in the element names
## may lead to names that are too long. Often, only the depth of a list
## element may only be important. To deal with this possibility a boolean
## option called `ForceNames` has been provided. ForceNames shall drop all
## original names in the lists and coerce simple numeric names which simply
## indicate the position of an element at the nth-level trunk as well as all
## preceding trunk numbers.
##
## Returns:
## LinearList: Named list.
##
LinearizeNestedList <- function(NList, LinearizeDataFrames=FALSE,
                                NameSep="/", ForceNames=FALSE) {
    ## Sanity checks:
    ##
    stopifnot(is.character(NameSep), length(NameSep) == 1)
    stopifnot(is.logical(LinearizeDataFrames), length(LinearizeDataFrames) == 1)
    stopifnot(is.logical(ForceNames), length(ForceNames) == 1)
    if (! is.list(NList)) return(NList)
    ##
    ## If no names on the top-level list coerce names. Recursion shall handle
    ## naming at all levels.
    ##
    if (is.null(names(NList)) | ForceNames == TRUE)
        names(NList) <- as.character(1:length(NList))
    ##
    ## If simply a dataframe deal promptly.
    ##
    if (is.data.frame(NList) & LinearizeDataFrames == FALSE)
        return(NList)
    if (is.data.frame(NList) & LinearizeDataFrames == TRUE)
        return(as.list(NList))
    ##
    ## Book-keeping code to employ a while loop.
    ##
    A <- 1
    B <- length(NList)
    ##
    ## We use a while loop to deal with the fact that the length of the nested
    ## list grows dynamically in the process of linearization.
    ##
    while (A <= B) {
        Element <- NList[[A]]
        EName <- names(NList)[A]
        if (is.list(Element)) {
            ##
            ## Before and After to keep track of the status of the top-level trunk
            ## below and above the current element.
            ##
            if (A == 1) {
                Before <- NULL
            } else {
                Before <- NList[1:(A - 1)]
            }
            if (A == B) {
                After <- NULL
            } else {
                After <- NList[(A + 1):B]
            }
            ##
            ## Treat dataframes specially.
            ##
            if (is.data.frame(Element)) {
                if (LinearizeDataFrames == TRUE) {
                    ##
                    ## `Jump` takes care of how much the list shall grow in this step.
                    ##
                    Jump <- length(Element)
                    NList[[A]] <- NULL
                    ##
                    ## Generate or coerce names as need be.
                    ##
                    if (is.null(names(Element)) | ForceNames == TRUE)
                        names(Element) <- as.character(1:length(Element))
                    ##
                    ## Just throw back as list since dataframes have no nesting.
                    ##
                    Element <- as.list(Element)
                    ##
                    ## Update names
                    ##
                    names(Element) <- paste(EName, names(Element), sep=NameSep)
                    ##
                    ## Plug the branch back into the top-level trunk.
                    ##
                    NList <- c(Before, Element, After)
                }
                Jump <- 1
            } else {
                NList[[A]] <- NULL
                ##
                ## Go recursive! :)
                ##
                if (is.null(names(Element)) | ForceNames == TRUE)
                    names(Element) <- as.character(1:length(Element))
                Element <- LinearizeNestedList(Element, LinearizeDataFrames,
                                               NameSep, ForceNames)
                names(Element) <- paste(EName, names(Element), sep=NameSep)
                Jump <- length(Element)
                NList <- c(Before, Element, After)
            }
        } else {
            Jump <- 1
        }
        ##
        ## Update book-keeping variables.
        ##
        A <- A + Jump
        B <- length(NList)
    }
    return(NList)
}



##' Display the SNPs with greatest marginal posterior probability of inclusion (MPPI)
##'
##' @title Best SNPs
##' @param d snpmod object or a list of snpmod objects
##' @param mppi.thr MPPI threshold, SNPs with MPPI>mppi.thr will be shown
##' @param pp.thr deprecated, alias for mppi.thr
##' @return subset of \code{snps(d)} data.frame including best SNPs
##' @author Chris Wallace
best.snps <- function(d,mppi.thr=pp.thr,pp.thr=0.1) {
  if(is.list(d))
    return(lapply(d,best.snps,mppi.thr=mppi.thr,pp.thr=pp.thr))
  if(!is(d,"snpmod"))
    stop("expected d to be a snpmod object, found ",class(d))
  tmp <- subset(d@snps, d@snps$Marg_Prob_Incl>mppi.thr)
  return(tmp[order(tmp$Marg_Prob_Incl,decreasing=TRUE),])
}
