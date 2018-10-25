#general diagnostic 

sig_est <- list(l4, b4)
res <- vector('list', length(sig_est))

for(sig in sig_est){
  #rep each curve
  modparams <- lapply(subsets, function(x) sub_genparams(listdf=x, est=sig))
  modresids <- lapply(modparams, function(x) lapply(x$fits, resid))
  modresids.unl <- unlist(modresids, recursive = FALSE)
  modresids.length <- lapply(modresids.unl, length)
  #finding rss total
  modrss <- lapply(modresids, function(x) lapply(x, function(x) sum(x^2)))
  #finding CT value
  mod1 <- lapply(subsets, function(x) modlist(x, model = l5))
  res1 <- lapply(mod1, function(x) getPar(x, type = "curve", cp = "cpD2", eff = "sliwin")) 
  res1CT <- lapply(res1, function(x) x[1,]) #first row is CT
  #specific CT +/- 2 cycles
  rssgrey <-list()
  for(i in 1:8){
    for(j in 1:12){
      indij = (12*(i-1))+j  
      uppercyc = floor(unlist(res1CT)[[indij]]+2)  #+2 cycles
      #klag-adjusted lower cycle
      lowercyc = ceiling(unlist(res1CT)[[indij]]-2) #-2cycles }
      #list of lists of rss grey region
      rssgrey[[indij]] <- sum(modresids.unl[[indij]][lowercyc:uppercyc]^2) 
    }
  }
  #autocorrelation and dw for amplification and residual values
  reg.amp <- vector('list', 8) ; reg.res <- vector('list', 8) #dynlm formula results
  dw.amp <- list() ; dw.res <- vector('list', 8) #dw-test values
  for(i in 1:8){
    reg.amp[[i]] <- lapply(subsets[[i]][ ,2:length(subsets[[i]])], 
                           function(y) dynlm(y ~ subsets[[i]][,1])) #all amp can run
    reg.res[[i]] <- lapply(modresids[[i]], function(y) try(dynlm(y ~ subsets[[i]][,1]), silent=TRUE))
    #residuals contain NA
    dw.amp[[i]] <- lapply(reg.amp[[i]], function(x) durbinWatsonTest(x)) #amplification values
  }
  for(i in 1:8){ #residuals
    for(j in 1:(length(subsets[[i]])-1)){ #NA resids cannot run
      dw.res[[i]][[j]] <- tryCatch({
        durbinWatsonTest(reg.res[[i]][[j]])
      }, error=function(e) {
        return(list(r=NA, dw=NA, p=NA))
      })
    } #replace errors with NA for r, dw, p
  }
  #extracting specific values
  dw.amp <- lapply(dw.amp, function(x) lapply(x, function(x) x[1:3]))
  dw.ampmat <- matrix(unlist(dw.amp), ncol = 3, byrow=TRUE)
  dw.res <- lapply(dw.res, function(x) lapply(x, function(x) x$r))
  dw.resmat <- matrix(unlist(dw.amp), ncol = 3, byrow=TRUE)
  #combining all into matrix
  empmat <- matrix(data=NA, nrow=96, ncol=10)
  empmat[,1] <- unique(unlist(lapply(subsets, colnames)))[-grep("Cycle", unique(unlist(lapply(subsets, colnames))))]
  empmat[,2] <- unlist(modrss)
  empmat[,3] <- unlist(res1CT)
  empmat[,4] <- unlist(rssgrey)
  empmat[,5:7] <- dw.ampmat
  empmat[,8:10] <- dw.resmat
  #dataframe output
  resstat <- data.frame(empmat)
  colnames(resstat) <- c("Replicate", "RSS", "CT", "RSSlocal", 
                         "amp.r", "amp.dw", "amp.p", 
                         "res.r", "res.dw", "res.p")
  
  if(identical(l4, sig)){
    res[[1]] <- resstat
  }
  if(identical(b4, sig)){
    res[[2]] <- resstat
  }
}
