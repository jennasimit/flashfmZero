#' @title Make table of fine-mapping results output from "flashfm" (flashfmZero or flashfm), "multiJAMd", or "JAMdwithGroups" 
#' @param FMout object output from "flashfm" (flashfmZero or flashfm related wrappers, like "FLASHFMZEROwithJAMd"), "multiJAMd", or "JAMdwithGroups" 
#' @param method text for fine-mapping method used, i.e. one of "flashfm" (flashfmZero or flashfm), "multiJAMd", or "JAMdwithGroups"
#' @param traitnames vector of trait names (same length as number of traits that were fine-mapped) and in same order as fine-mapped traits
#' @param cred level used for credible set construction; default is 0.99
#' @return Table where each row corresponds to one trait, listing credible set size, details of variant with maximum MPP, variants with MPP>0.90 and snp groups coinciding with the variant(s)
#' @export
#' @author Jenn Asimit
FMsummary_table <- function(FMout,method="flashfm",traitnames=NULL,cred=0.99) {
# method ="flashmfm" for output from flashfm; method="multiJAMd"; method="JAMdwithGroups"
 
if(method=="flashfm") {
 FMoutCS <- allcredsetsPP(FMout$mpp.pp,cred=cred)
# traitnames <- names(FMoutCS$fm)
 CSsizes1 <-sapply(FMoutCS$fm, nrow) # "1" for single-trait
 CSsizesM <-sapply(FMoutCS$flashfm, nrow) # "M" for multi-trait
 topSNP1 <- sapply(FMoutCS$fm, function(x) x[1,"SNP"]) # snps are sorted by MPP, so first snp had max MPP - output snp name and its MPP
 topSNPM <- sapply(FMoutCS$flashfm, function(x) x[1,"SNP"])
 topMPP1 <- sapply(FMoutCS$fm, function(x) x[1,"MPP"]) # snps are sorted by MPP, so first snp had max MPP - output snp name and its MPP
 topMPPM <- sapply(FMoutCS$flashfm, function(x) x[1,"MPP"])

 snpGroups <-  FMout$snpGroups
 group1 <- groupIDs.fn(snpGroups$groups.fm,topSNP1) 
 groupM <- groupIDs.fn(snpGroups$groups.flashfm,topSNPM)
 group1PP <- sapply(1:length(group1),function(x) FMout$mpp.pp$MPPg[[x]][group1[x],1]) 
 groupMPP <- sapply(1:length(groupM),function(x) FMout$mpp.pp$MPPg[[x]][groupM[x],2])
 group1size <- as.vector(snpGroups$group.sizes[1,group1])
 group1size <- as.vector(unname(group1size),mode="numeric")
 groupMsize <- snpGroups$group.sizes[2,groupM]
 groupMsize <- as.vector(unname(groupMsize),mode="numeric")
 
 SNP1.gt90 <- sapply(FMout$mpp.pp$MPP, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- rownames(x)[ind]; out <- paste(out,collapse=",")};return(out) } )
 SNPM.gt90 <- sapply(FMout$mpp.pp$MPP, function(x) { ind <- which(x[,2]>0.90); out <- "NA"; if(length(ind)>0) {out <- rownames(x)[ind]; out <- paste(out,collapse=",")};return(out) } )
 MPP1.gt90 <- sapply(FMout$mpp.pp$MPP, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- x[ind,1]; out <- paste(out,collapse=",")};return(out) } )
 MPPM.gt90 <- sapply(FMout$mpp.pp$MPP, function(x) { ind <- which(x[,2]>0.90); out <- "NA"; if(length(ind)>0) {out <- x[ind,2]; out <- paste(out,collapse=",")};return(out) } )
 
 group1.gt90 <- sapply(FMout$mpp.pp$MPPg, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- rownames(x)[ind]; out <- paste(out,collapse=",")};return(out) } )
 groupM.gt90 <- sapply(FMout$mpp.pp$MPPg, function(x) { ind <- which(x[,2]>0.90); out <- "NA"; if(length(ind)>0) {out <- rownames(x)[ind]; out <- paste(out,collapse=",")};return(out) } )
 MPPg1.gt90 <- sapply(FMout$mpp.pp$MPPg, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- x[ind,1]; out <- paste(out,collapse=",")};return(out) } )
 MPPgM.gt90 <- sapply(FMout$mpp.pp$MPPg, function(x) { ind <- which(x[,2]>0.90); out <- "NA"; if(length(ind)>0) {out <- x[ind,2]; out <- paste(out,collapse=",")};return(out) } )

 group1size.gt90 <- sapply(group1.gt90,function(x) {out <- "NA"; if(!is.na(x)){ g <- unlist(strsplit(x,",")); l <- sapply(g,function(x) length(snpGroups$groups.fm[[x]])); out <- paste(l,collapse=",")}; return(out)})
 groupMsize.gt90 <- sapply(groupM.gt90,function(x) {out <- "NA"; if(!is.na(x)){ g <- unlist(strsplit(x,",")); l <- sapply(g,function(x) length(snpGroups$groups.flashfm[[x]])); out <- paste(l,collapse=",")}; return(out)})
 
 out1 <- data.frame(Approach="JAM-latent-factor",Trait=traitnames,CS99size=CSsizes1,topSNP=topSNP1,topMPP=topMPP1,topGroup=group1,topGroupSize=group1size,topGroupMPP=group1PP, SNP_gt90=SNP1.gt90, MPP_gt90=MPP1.gt90, group_gt90=group1.gt90, group_gt90_size=group1size.gt90, MPPg_gt90= MPPg1.gt90)
 out1 <- out1[order(out1$Trait),]
 outM <- data.frame(Approach="flashfm-latent-factor",Trait=traitnames,CS99size=CSsizesM,topSNP=topSNPM,topMPP=topMPPM,topGroup=groupM,topGroupSize=groupMsize,topGroupMPP=groupMPP, SNP_gt90=SNPM.gt90, MPP_gt90=MPPM.gt90,group_gt90=groupM.gt90, group_gt90_size=groupMsize.gt90, MPPg_gt90= MPPgM.gt90)
 outM <- outM[order(outM$Trait),]
 out <- rbind(out1,outM)
 rownames(out) <- NULL
 out$topGroup <- paste0(out$topGroup,"_lat")
 
}

if(method=="multiJAMd"){
 FMoutCS <- multiJAMdCS(FMout, cred=cred)
# traitnames <- names(FMoutCS$fm)
 CSsizes1 <-sapply(FMoutCS$fm, nrow) 
 topSNP1 <- sapply(FMoutCS$fm, function(x) x[1,"SNP"]) # snps are sorted by MPP, so first snp had max MPP - output snp name and its MPP
 topMPP1 <- sapply(FMoutCS$fm, function(x) x[1,"MPP"]) # snps are sorted by MPP, so first snp had max MPP - output snp name and its MPP

 snpGroups <-  FMout$snpGroups
 group1 <- groupIDs.fn(snpGroups,topSNP1) 
 group1PP <- sapply(1:length(group1),function(x) FMout$mpp.pp$MPPg[[x]][group1[x],"MPP"]) 
 group1size <- sapply(group1,function(x) length(snpGroups[[x]]))
 group1size <- as.vector(unname(group1size),mode="numeric")

 SNP1.gt90 <- sapply(FMout$mpp.pp$MPP, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- rownames(x)[ind]; out <- paste(out,collapse=",")};return(out) } )
 MPP1.gt90 <- sapply(FMout$mpp.pp$MPP, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- x[ind,1]; out <- paste(out,collapse=",")};return(out) } )

 group1.gt90 <- sapply(FMout$mpp.pp$MPPg, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- rownames(x)[ind]; out <- paste(out,collapse=",")};return(out) } )
 MPPg1.gt90 <- sapply(FMout$mpp.pp$MPPg, function(x) { ind <- which(x[,1]>0.90); out <- "NA"; if(length(ind)>0) {out <- x[ind,1]; out <- paste(out,collapse=",")};return(out) } )

 group1size.gt90 <- sapply(group1.gt90,function(x) {out <- "NA"; if(!is.na(x)){ g <- unlist(strsplit(x,",")); l <- sapply(g,function(x) length(snpGroups[[x]])); out <- paste(l,collapse=",")}; return(out)})

 out <- data.frame(Approach="JAM-observed-trait",Trait=traitnames,CS99size=CSsizes1,topSNP=topSNP1,topMPP=topMPP1,topGroup=group1,topGroupSize=group1size,topGroupMPP=group1PP,SNP_gt90=SNP1.gt90, MPP_gt90=MPP1.gt90, group_gt90=group1.gt90, group_gt90_size=group1size.gt90,MPPg_gt90= MPPg1.gt90)
 out <- out[order(out$Trait),]
 rownames(out) <- NULL
  out$topGroup <- paste0(out$topGroup,"_obs")
}

if(method=="JAMdwithGroups") {
 FMoutCS <- FMout$CS
 CSsizes1 <- nrow(FMoutCS)
 topSNP1 <- FMoutCS$snp[1]
 topMPP1 <- FMoutCS$MPP[1]
 
 group1 <- FMoutCS$group[1]
 group1PP <- FMout$mpp.pp$MPPg[group1,"MPP"]
 group1size <- length(FMout$snpGroups[[group1]])
 
 ind.gt90 <- which(FMout$mpp.pp$MPP>0.90)
 SNP1.gt90 <- MPP1.gt90 <- "NA"
 if(length(ind.gt90)>0) {
  SNP1.gt90 <- paste(rownames(FMout$mpp.pp$MPP)[ind.gt90],collapse=",")
  MPP1.gt90 <- paste(FMout$mpp.pp$MPP[ind.gt90,1],collapse=",")
 }
 
 ind.gt90 <- which(FMout$mpp.pp$MPPg>0.90)
 group1.gt90 <- MPPg1.gt90 <- group1size.gt90 <- "NA"
 if(length(ind.gt90)>0) {
  group1.gt90 <- paste(rownames(FMout$mpp.pp$MPPg)[ind.gt90],collapse=",")
  MPPg1.gt90 <- paste(FMout$mpp.pp$MPPg[ind.gt90,1],collapse=",")
 snpGroups <- FMout$snpGroups
 group1size.gt90 <- sapply(group1.gt90,function(x) {out <- "NA"; if(!is.na(x)){ g <- unlist(strsplit(x,",")); l <- sapply(g,function(x) length(snpGroups[[x]])); out <- paste(l,collapse=",")}; return(out)})
}
 out <- data.frame(Approach="JAM-latent-factor",Trait=traitnames,CS99size=CSsizes1,topSNP=topSNP1,topMPP=topMPP1,topGroup=group1,topGroupSize=group1size,topGroupMPP=group1PP,SNP_gt90=SNP1.gt90, MPP_gt90=MPP1.gt90, group_gt90=group1.gt90, group_gt90_size=group1size.gt90,MPPg_gt90= MPPg1.gt90)
 out <- out[order(out$Trait),]
 rownames(out) <- NULL
  out$topGroup <- paste0(out$topGroup,"_lat")
 
 
}


return(out)
}


#' @title Make table of fine-mapping results output for latent factors and observed traits (wrapper for FMsummary_table)
#' @param FMobs observed trait fine-mapping object output from "multiJAMd"
#' @param FMlatent latent factor object output from "flashfm" (flashfmZero or flashfm related wrappers, like "FLASHFMZEROwithJAMd") or "JAMdwithGroups" 
#' @param fm_traits_ob vector of observed trait names (same length as number of traits that were fine-mapped) and in same order as fine-mapped traits
#' @param fm_traits_latent vector of latent factor names (same length as number of factors that were fine-mapped) and in same order as fine-mapped factors
#' @param cred level used for credible set construction; default is 0.99
#' @param regions data.frame with region details, where the columns are in order: "chromosom"e, "start"", "end"
#' @param array_index row number of regions data.frame 
#' @return Table where each row corresponds to one trait, listing credible set size, details of variant with maximum MPP, variants with MPP>0.90 and snp groups coinciding with the variant(s)
#' @export
#' @author Jenn Asimit
FMsummary_table_general <- function(FM_obs, FM_latent, fm_traits_ob, fm_traits_latent,cred=0.99, regions=NULL,array_index=NULL){
 out1 <- FMsummary_table(FM_obs,method="multiJAMd",traitnames=fm_traits_obs,cred=cred)
 if(names(FM_latent)[1] == "mpp.pp") {
  outM <- FMsummary_table(FM_latent,method="flashfm",traitnames=fm_traits_latent,cred=cred)
 } else { outM <- FMsummary_table(FM_latent,method="JAMdwithGroups",traitnames=fm_traits_latent,cred=cred) }
 out <- rbind(out1,outM)
 if(!is.null(regions) & !is.null(array_index)) { out <- cbind(regions[rep(array_index,nrow(out)),],out)} 
rownames(out) <- NULL
return(out)
}


# Example use:
#out <- NULL
#regions <- read.table("region_list.tsv",header=TRUE)
#for(array_index in 1:nrow(regions)) {
#  if(file.exists(paste0(array_index,".Rdata"))){
#  load(paste0(array_index,".Rdata"))
#  out1 <- FMsummary_table_general(FM_obs, FM_latent, fm_traits_ob, fm_traits_latent,regions=regions,array_index=array_index,cred=0.99) 
#  out <- rbind(out,out1)
#  }
#}


