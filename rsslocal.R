allresids <- function(subs, params.sig, params.lstar, sig=l5){
  ##extracting results from matrix
  res.sig <- gen_results(subs, params.sig)
  res.lstar <- gen_results(subs, params.lstar)
  sigCT <- res.sig$CT
  lstarCT <- res.lstar$CT
  #branching on its own
  branch.res <- branchfunc(F)
  branchCT <- branch.res$CT
  #branchCT.mat <- matrix(branch.res$CT, nrow=8, ncol=12, byrow=TRUE)
  #branchCT <- list()
  #for(i in 1:8){
  #  branchCT[[i]] <- unlist(branchCT.mat[i,])
  #}
  ##finding resids 
  #resids with NA get filtered out
  residsNA <- function(x){
    tryCatch({ resid(x) }, error = function(e){NA})
  }
  #sigmoidal
  modparams.sig <- lapply(subs, function(x) sub_genparams(listdf=x, est=sig))
  modresids.sig <- lapply(modparams.sig, function(x) lapply(x$fits, residsNA)) #for non-models, only one NA
  #lstar
  modresids.lstar <- lapply(params.lstar, function(x) lapply(x$fits, residsNA)) #for non-models, only one NA
  #branching
  modresids.branch <- branchCTunl$residuals
  #unlisting residuals
  oneNA <- function(modresids){
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        if(all(is.na(modresids[[i]][[j]])) == TRUE){
          modresids[[i]][[j]] <- rep(NA, cyclength[[i]])
        }
        else{}
      }
    }
    return(modresids)
  }
  modresids.sigunl <- oneNA(modresids.sig) 
  modresids.lstarunl <- oneNA(modresids.lstar) 
  modresids.branchunl <- oneNA(modresids.branch)
  ##unlisting resids
  modresids.sigunl <- unlist(modresids.sigunl, recursive = FALSE)
  modresids.lstarunl <- unlist(modresids.lstarunl, recursive = FALSE)
  modresids.branchunl <- unlist(modresids.branchunl, recursive = FALSE)
  
  ####+/- 2 cycles SIGMOIDAL
  sigrsslocal <- function(res, modresids.unl){
    rssgrey.sig <-list()
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        indij = ((replength/sublength)*(i-1))+j  
        uppercyc = floor(as.numeric(as.character(res[[indij]]))+2)  #+2 cycles
        #klag-adjusted lower cycle
        lowercyc = ceiling(as.numeric(as.character(res[[indij]]))-2) #-2cycles }
        #list of lists of rss grey region
        rssgrey.sig[[indij]] <- tryCatch({ sum(modresids.unl[[indij]][lowercyc:uppercyc]^2) },
                                         error = function(e) { NA })
      }
    }
    return(unlist(rssgrey.sig))
  }
  
  ####+/- 2 cycles LSTAR
  
  findlstarres <- function(thr, params){
    mod2unl <- unlist(unlist(thr, recursive=FALSE))
    #add NA to those that are outside lower bound first
    mod2unl[which(mod2unl < (0.5+klag*mdim+2))] <- NA
    #rss for red box
    rsslstarr <- vector('list', sublength) 
    for(i in 1:(sublength)){
      for(j in 1:(replength/sublength)){
        #combined indicator with i and j
        indij = ((replength/sublength)*(i-1))+j  
        #all NA values
        if(sum(is.na(mod2unl[[indij]])) > 0){
          #this does not go through, will have to run this at the end
          rsslstarr[[i]][[j]] <- NA
        }
        #filter out low values, change to NA unlistcycCTthr
        else if(sum(isTRUE(mod2unl[[indij]] > (cyclength[[i]]-(klag*mdim+1)))) > 0){ #strict upper bound
          uppercyc = cyclength[[i]]-(klag*mdim)
          lowercyc = ceiling(as.numeric(as.character(mod2unl[[indij]]))-(klag*mdim)-2)}
        else{    #(sum(is.na(unlistcycCTthr[[indij]])) == 0){ 
          uppercyc = floor(as.numeric(as.character(mod2unl[[indij]]))-(klag*mdim)+2)  #+2 cycles
          #klag-adjusted lower cycle
          lowercyc = ceiling(as.numeric(as.character(mod2unl[[indij]]))-(klag*mdim)-2)} #-2cycles 
        #list of lists of rss red region
        rsslstarr[[i]][[j]] <- sum(params[[i]]$fits[[j]]$residuals[lowercyc:uppercyc]^2)
        #rsslstarr provies RSS for NA values, we undo this here (rewrite above, need both)
        if(sum(is.na(mod2unl[[indij]])) > 0){ 
          rsslstarr[[i]][[j]] <- NA }
        else{}
      }
    }
    return(unlist(rsslstarr))
  }
  #rss local for sig.sig
  rsslocal.sig <- sigrsslocal(sigCT, modresids.sigunl)
  rsslocal.sig.lstar <- sigrsslocal(lstarCT, modresids.sigunl)
  rsslocal.sig.branch <- sigrsslocal(branchCT, modresids.sigunl)
  rsslocal.sig.all <- matrix(c(rsslocal.sig, rsslocal.sig.lstar, rsslocal.sig.branch), ncol=3)
  #unlisted 1:(replength) of rss grey region
  rsslocal.lstar <- findlstarres(lstarCT, lstarres.fits)
  rsslocal.lstar.sig <- findlstarres(sigCT, lstarres.fits)
  rsslocal.lstar.branch <- findlstarres(branchCT, lstarres.fits)
  rsslocal.lstar.all <- matrix(c(rsslocal.lstar, rsslocal.lstar.sig, rsslocal.lstar.branch), ncol=3)
  ##resids for branching
  rsslocal.branch <- sigrsslocal(branchCT, modresids.branchunl)
  rsslocal.branch.sig <- sigrsslocal(sigCT, modresids.branchunl)
  rsslocal.branch.lstar <- sigrsslocal(lstarCT, modresids.branchunl)
  rsslocal.branch.all <- matrix(c(rsslocal.branch, rsslocal.branch.sig, rsslocal.branch.lstar), ncol=3)
  
  return(list(sig = rsslocal.sig.all, lstar = rsslocal.lstar.all, branch = rsslocal.branch.all))
}

meplz <- allresids(subsets, modparams, lstarres.fits)



#finding CT values
gap.sig <- gen_results(subsets, modparams)
gap.lstar <- gen_results(subsets, lstarres.fits)
sigCT <- gap.sig$CT
lstarCT <- gap.lstar$CT

#binomial branching CT
F = t(do.call(cbind, lapply(subsets,function(x) x[-1])))
branchCTunl <- branchfunc(F)
branchCTunl.mat <- matrix(branchCTunl$CT, nrow=8, ncol=12, byrow=TRUE)
branchCT <- list()
for(i in 1:8){
  branchCT[[i]] <- unlist(branchCTunl.mat[i,])
}

##finding resids 
#resids with NA get filtered out
residsNA <- function(x){
  tryCatch({ resid(x) }, error = function(e){NA})
}
#sigmoidal
modparams.sig <- lapply(subsets, function(x) sub_genparams(listdf=x, est=l5))
modresids.sig <- lapply(modparams.sig, function(x) lapply(x$fits, residsNA)) #for non-models, only one NA
#lstar
modresids.lstar <- lapply(lstarres.fits, function(x) lapply(x$fits, residsNA)) #for non-models, only one NA
#branching
modresids.branch <- branchCTunl$residuals

#unlisting residuals
oneNA <- function(modresids){
  for(i in 1:sublength){
    for(j in 1:(replength/sublength)){
      if(all(is.na(modresids[[i]][[j]])) == TRUE){
        modresids[[i]][[j]] <- rep(NA, cyclength[[i]])
      }
      else{}
    }
  }
  return(modresids)
}
modresids.sigunl <- oneNA(modresids.sig) 
modresids.lstarunl <- oneNA(modresids.lstar) 
modresids.branchunl <- oneNA(modresids.branch)
##unlisting resids
modresids.sigunl <- unlist(modresids.sigunl, recursive = FALSE)
modresids.lstarunl <- unlist(modresids.lstarunl, recursive = FALSE)
modresids.branchunl <- unlist(modresids.branchunl, recursive = FALSE)
