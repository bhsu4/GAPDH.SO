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



