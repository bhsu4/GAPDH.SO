res_qc <- function(getfiles, orgdata, est, res_out){

  files <- getfiles
  targnames <- unique(files) #all targets
  targnames <- gsub(".Rda", "", targnames) #remove .Rda to match later
  sampnames <- sort(unique(c(orgdata$SampleID)))
  #repeat all sample names for the chosen targets in files
  tmp <- rep(sampnames, length(files))
  #creating result data.frame
  targlength <- length(targnames) #length of all targets
  replength <- length(unique(orgdata$SampleID)) #sum(lengths(tst)) #length of total reps for each target
  res <- data.frame(
    #target categories
    TargetName = rep(targnames, each = replength), SampleID = rep(sampnames, targlength), 
    Group = gsub("_." , "", tmp),
    #r-squared quality scores
    Rsq = rep(NA, targlength * replength)
  )
  #start of calculating r-squared quality scores
  foreach(k = 1:length(files)) %do% {    #(k in 1:length(files)){  #foreach (k = 1:length(files)) %do% {
    load(file = files[[k]])
    print(files[[k]])
    try <- unlist.genparams(tst)
    #finding quality score
    res_pcr <- lapply(try, function(x) sub_genparams(est= est, listdf = x))
    #tryCatch non-sigmoidal fits to return NA
    effr2 <- function(x) tryCatch({efficiency(x, plot=FALSE)$Rsq}, error = function(e) NA)  
    res_r2 <- lapply(res_pcr, function(x) lapply(x$fits, function(x) effr2(x)))
    #indexing 
    ind2 <- length(unique(orgdata$SampleID))*k  ; ind1 <- ind2-(length(unique(orgdata$SampleID))-1)
    res[ind1:ind2, "Rsq"] <- unlist(res_r2)
  }
  #left join outputs to create qc column with matched IDs
  res_out[,"id"] <- paste(res_out$TargetName, res_out$SampleID, sep="-")
  res[,"id"] <- paste(res$TargetName, res$SampleID, sep="-")
  res_merge <- merge(res_out, res, by="id", all.x=TRUE)
  return(res_merge)
}

gen_ctqc <- function(res){
  #parameter lengths
  replength <- length(unique(res$SampleID.x))
  grouplength <- length(unique(res$Group.x))
  #empty matrix for ct and rsq 
  res_ct <- matrix(NA, length(unique(res$TargetName.x)), replength)
  res_rsq <- matrix(NA, length(unique(res$TargetName.x)), replength)
  #creating ct and rsq matrix
  for(k in 1:length(unique(res$TargetName.x))){
    ind2 <- replength*k  ; ind1 <- ind2-(replength-1)
    res_curr <- res[ind1:ind2, c("SampleID.x", "ct", "Rsq")]
    sampids <- res[ind1:ind2, c("SampleID.x")] #all unique sample IDs (KW)
    sampids_ord <- order(match(sampids, paste0(mixedsort(gsub( "_.*$", "", 
                         res[ind1:ind2, c("SampleID.x")])), "_", 1:(replength/grouplength))))
    ct_val <- t(res_curr[sampids_ord,]$ct) #correct order, mixed sorting
    rsq_val <- t(res_curr[sampids_ord,]$Rsq)
    #creating matrix for ct
    res_ct[k,] <- ct_val
    res_rsq[k,] <- rsq_val
    #matrix names
    rownames(res_ct) <- unique(res$TargetName.x)
    rownames(res_rsq) <- unique(res$TargetName.x)
    colnames(res_ct) <- as.vector(paste0(rep(paste0("KW", 1:grouplength), 
                                  each = (replength/grouplength)), ":", 1:(replength/grouplength))) #KW(n)
    colnames(res_rsq) <- as.vector(paste0(rep(paste0("KW", 1:grouplength), 
                                   each = (replength/grouplength)), ":", 1:(replength/grouplength)))
  }
  res_qpcr <- list(ct=res_ct, qc=res_rsq)
  return(res_qpcr)
}






