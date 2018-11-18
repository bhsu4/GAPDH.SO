branchfunc <- function(subs, Fer, plot=FALSE){
  #setting lengths
  sublength = length(subs)
  repnames <- unique(unlist(lapply(subs, names)))[!grepl("Cycle", unique(unlist(lapply(subs, names))))]
  replength = length(repnames)
  cyclength <- unlist(lapply(subs, nrow))
  #start
  F = Fer
  a1 <- dim(F)
  r <- a1[1]
  n <- a1[2]
  tau1 <- 1.02
  tau2 <- 1.02
  ctau1 <- matrix(0, r, 1)
  ctau2 <- matrix(0, r, 1)
  for (i in 1:r) {
    if(F[i, 1] == 0){
      ctau1[i] <- min(which(F[i, 3:n]/F[i, 2:(n - 1)] >= tau1))
      if (any(F[i, (ctau1[i] + 1):n]/F[i, ctau1[i]:(n - 1)] <= tau2)) {
        ctau2[i] <- min(which(F[i, (ctau1[i] + 1):n]/F[i, ctau1[i]:(n - 1)] <= tau2)) + ctau1[i]
      }
      else {
        ctau2[i] = n
      }
    }
    ctau1[i] <- min(which(F[i, 2:n]/F[i, 1:(n - 1)] >= tau1))
    if (any(F[i, (ctau1[i] + 1):n]/F[i, ctau1[i]:(n - 1)] <= tau2)) {
      ctau2[i] <- min(which(F[i, (ctau1[i] + 1):n]/F[i, ctau1[i]:(n - 1)] <= tau2)) + ctau1[i]
    }
    else {
      ctau2[i] = n
    }
  }
  #}
  yn <- matrix(0, r, 1)
  yn_1 <- matrix(0, r, 1)
  p_tilde <- matrix(0, r, 1)
  mu_b <- matrix(0, r, 1)
  var_mu_b <- matrix(0, r, 1)
  for (i in 1:r) {
    yn[i] <- sum(F[i, ctau1[i]:ctau2[i]])
    yn_1[i] <- sum(F[i, ctau1[i]:(ctau2[i] - 1)])
    p_tilde[i] <- (yn[i] - F[i, ctau1[i]] - yn_1[i])/yn_1[i]
    mu_b[i] <- p_tilde[i]/((1 + p_tilde[i])^n) * yn[i] * (1 - (1 + p_tilde[i])^(ctau1[i] - n))^(-1)
  }
  
  n_an <- matrix(0, r, n) ; cdivas <- matrix(0, r, 1) ; appF <- matrix(0, r, n)
  n_an[,1] <- mean(mu_b)
  for(i in 1:r){
    for(j in 2:n){
      n_an[i,j] <- n_an[i, j-1] + rbinom(1, round(n_an[i, j-1]), p_tilde[i])
    }
    cdivas[i,] = (n_an[i, round((ctau1[i,]+ctau2[i,])/2)]/F[i,round((ctau1[i,]+ctau2[i,])/2)])/(9.1*10^11)
    constant = cdivas * (9.1*10^11)
    appF[i,] = n_an[i,]/constant[i,]
  }
  
  backconst <- matrix(0, r, n) ; res_aq <- matrix(0, r, n)
  for(i in 1:r){
    backconst[i,] = cumsum(F[i,])/F[i,]
    res_aq[i,] <- cumsum(appF[i,])/backconst[i,]
  }
  if(plot){  
    for(i in 1:r){
      plot(x=1:40, y=res_aq[i,], col=2)
      points(x=1:40, y=F[i,])
    }
  }
  resid.res <- vector('list', sublength)
  for(i in 1:sublength){
    for(j in 1:(replength/sublength)){
      indij = (replength/sublength)*(i-1)+j
      resid.res[[i]][[j]] <- res_aq[indij,] - F[indij,]
    }
  }
  branchCT <- list()
  for(i in 1:r){
    branchCT[i] <- quantile(ctau1[i,]:ctau2[i,])[2] #1st quartile
  }
  return(list(CT = branchCT, dat = res_aq, residuals = resid.res, exp = cbind(ctau1, ctau2)))
}

allresids <- function(subs, params.sig, params.lstar, Fer, sig=l5){
  #lengths
  sublength = length(subs)
  repnames <- unique(unlist(lapply(subs, names)))[!grepl("Cycle", unique(unlist(lapply(subs, names))))]
  replength = length(repnames)
  cyclength <- unlist(lapply(subs, nrow))
  ##extracting results from matrix
  res.sig <- gen_results(subs, params.sig)
  res.lstar <- gen_results(subs, params.lstar)
  sigCT <- res.sig$CT
  lstarCT <- res.lstar$CT
  #branching on its own
  branch.res <- branchfunc(subs, Fer)
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
  modresids.branch <- branch.res$residuals
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
    rssgrey.sig <-list() ; rssgrey.sigm <- list()
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        indij = ((replength/sublength)*(i-1))+j  
        uppercyc = floor(as.numeric(as.character(res[[indij]]))+2)  #+2 cycles
        #klag-adjusted lower cycle
        lowercyc = ceiling(as.numeric(as.character(res[[indij]]))-2) #-2cycles }
        #list of lists of rss grey region
        rssgrey.sig[[indij]] <- tryCatch({ sum(modresids.unl[[indij]][lowercyc:uppercyc]^2) },
                                         error = function(e) { NA })
        rssgrey.sigm[[indij]] <- rssgrey.sig[[indij]]/(uppercyc-lowercyc+1)
        
      }
    }
    return(matrix(c(unlist(rssgrey.sig), unlist(rssgrey.sigm)), ncol=2))
  }
  
  ####+/- 2 cycles LSTAR
  
  findlstarres <- function(thr, params){
    mod2unl <- unlist(unlist(thr, recursive=FALSE))
    #add NA to those that are outside lower bound first
    mod2unl[which(mod2unl < (0.5+klag*mdim+2))] <- NA
    #rss for red box
    rsslstarr <- lapply(vector('list', sublength), function(x) rep(0,(replength/sublength)))
    rsslstarrm <- lapply(vector('list', sublength), function(x) rep(0,(replength/sublength)))
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        #combined indicator with i and j
        indij = ((replength/sublength)*(i-1))+j  
        #all NA values
        if(sum(is.na(mod2unl[[indij]])) > 0){
          #this does not go through, will have to run this at the end
          rsslstarr[[i]][[j]] <- NA
        }
        #filter out low values, change to NA unlistcycCTthr
        else if(isTRUE(mod2unl[[indij]] > (cyclength[[i]]-(klag*mdim+1))) ==TRUE){ #strict upper bound
          uppercyc = cyclength[[i]]-(klag*mdim)
          lowercyc = ceiling(as.numeric(as.character(mod2unl[[indij]]))-(klag*mdim)-2)
          #list of lists of rss red region
          rsslstarr[[i]][[j]] <- sum(params[[i]]$fits[[j]]$residuals[lowercyc:uppercyc]^2)}
        else{    #(sum(is.na(unlistcycCTthr[[indij]])) == 0){ 
          uppercyc = floor(as.numeric(as.character(mod2unl[[indij]]))-(klag*mdim)+2)  #+2 cycles
          #klag-adjusted lower cycle
          lowercyc = ceiling(as.numeric(as.character(mod2unl[[indij]]))-(klag*mdim)-2)
          #list of lists of rss red region
          rsslstarr[[i]][[j]] <- sum(params[[i]]$fits[[j]]$residuals[lowercyc:uppercyc]^2)} #-2cycles 
          rsslstarrm[[i]][[j]] <- rsslstarr[[i]][[j]]/(uppercyc-lowercyc+1)
      }
    }
    return(matrix(c(unlist(rsslstarr), unlist(rsslstarrm)), ncol=2))
  }
  
  ####+/- exp phase BRANCHING
  branchlocal <- function(branch.res){
    modresids.branch <- branch.res$residuals
    modresids.branchunl <- oneNA(modresids.branch)
    modresids.branchunl <- unlist(modresids.branchunl, recursive = FALSE)
    rssgrey.branch <-list() ; rssgrey.branchm <- list()
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        indij = ((replength/sublength)*(i-1))+j  
        uppercyc = ctau1[indij]  #+2 cycles
        #klag-adjusted lower cycle
        lowercyc = ctau2[indij] #-2cycles }
        #list of lists of rss grey region
        rssgrey.branch[[indij]] <- tryCatch({ sum(modresids.branchunl[[indij]][lowercyc:uppercyc]^2) },
                                         error = function(e) { NA })
        rssgrey.branchm[[indij]] <- rssgrey.branch[[indij]]/(uppercyc-lowercyc+1)
      }
    }
    return(matrix(c(unlist(rssgrey.sig), unlist(rssgrey.sigm)), ncol=2))
  }
  
  
  #rss local for sig.sig
  rsslocal.sig <- sigrsslocal(sigCT, modresids.sigunl)
  #unlisted 1:(replength) of rss grey region
  rsslocal.lstar <- findlstarres(lstarCT, lstarres.fits)
  ##resids for branching
  rsslocal.branch <- sigrsslocal(branchCT, modresids.branchunl)
  #all together
  rss.locals <- data.frame(matrix(c(rsslocal.sig, rsslocal.lstar, rsslocal.branch), ncol=6))
  colnames(rss.locals) <- c("RSS_Sig", "RSSLocal_Sig", "RSS_LSTAR", "RSSLocal_LSTAR", "RSS_Branch", "RSSLocal_Branch")
  
  return(rss.locals)
}


#GAPDH.SO TESTING
#outside params generated
subsetslog = lapply(subsets, function(x) log10(x))
for(i in 1:8){
  subsetslog[[i]][,1] <- 1:40
}
modparams <- lapply(subsets, function(x) sub_genparams(listdf=x, est=l5))
lstarres <- vector("list", 8) #empty list of subs: ABCD,..etc
cyclength <- unlist(lapply(subsets, nrow))
for(i in 1:8){ #running LSTAR model
  for(j in 2:13){
    #results for LSTAR model 
    lstarres[[i]][[j-1]] <- tryCatch({
      tsDyn::lstar(subsets[[i]][,j], m=mdim, d=klag)}, #d = lag found through AIC
      error=function(e) list(fitted.values=rep(NA, (cyclength[[i]]-(klag*mdim))), 
                             residuals=rep(NA, (cyclength[[i]]-(klag*mdim))),
                             model.specific=list(coefficients = rep(NA, (4+2*mdim)))))
    #if error, output NA, (no fit) w/ fitted values and residual rep length times minus klag
  }
}
lstarres.fits <- lapply(lstarres, function(x) list(fits=x))

F = t(do.call(cbind, lapply(subsets,function(x) x[-1])))
meplz <- allresids(subsets, modparams, lstarres.fits, Fer = F)

do.call(cbind, meplz)

#try bad data set mircomp
#not nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst)
#outside params generated
modparams <- lapply(try, function(x) sub_genparams(listdf=x, est=l5))
lstarres <- vector("list", 10) #empty list of subs: ABCD,..etc
cyclength <- unlist(lapply(try, nrow))
for(i in 1:10){ #running LSTAR model
  for(j in 2:5){
    #results for LSTAR model 
    lstarres[[i]][[j-1]] <- tryCatch({
      tsDyn::lstar(try[[i]][,j], m=mdim, d=klag)}, #d = lag found through AIC
      error=function(e) list(fitted.values=rep(NA, (cyclength[[i]]-(klag*mdim))), 
                             residuals=rep(NA, (cyclength[[i]]-(klag*mdim))),
                             model.specific=list(coefficients = rep(NA, (4+2*mdim)))))
    #if error, output NA, (no fit) w/ fitted values and residual rep length times minus klag
  }
}
lstarres.fits <- lapply(lstarres, function(x) list(fits=x))
F = t(do.call(cbind, lapply(try, function(x) x[-1][1:40,])))
meplz2 <- allresids(try, modparams, lstarres.fits, Fer = F)


#nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams.Rn(tst)
#nice
modparams.good <- lapply(try.good, function(x) sub_genparams(listdf=x, est=l5))
lstarres.good <- vector("list", 10) #empty list of subs: ABCD,..etc
cyclength <- unlist(lapply(try.good, nrow))
for(i in 1:10){ #running LSTAR model
  for(j in 2:5){
    #results for LSTAR model 
    lstarres.good[[i]][[j-1]] <- tryCatch({
      tsDyn::lstar(try.good[[i]][,j], m=mdim, d=klag)}, #d = lag found through AIC
      error=function(e) list(fitted.values=rep(NA, (cyclength[[i]]-(klag*mdim))), 
                             residuals=rep(NA, (cyclength[[i]]-(klag*mdim))),
                             model.specific=list(coefficients = rep(NA, (4+2*mdim)))))
    #if error, output NA, (no fit) w/ fitted values and residual rep length times minus klag
  }
}
lstarres.fits <- lapply(lstarres.good, function(x) list(fits=x))
F = t(do.call(cbind, lapply(try.good, function(x) x[-1][1:40,])))
meplz3 <- allresids(try.good, modparams.good, lstarres.fits, Fer = F)


#testing on another nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-520e_001119.Rda")
try.nice <- unlist.genparams.Rn(tst)
modparams <- lapply(try.nice, function(x) sub_genparams(listdf=x, est=l5))
lstarres <- vector("list", 10) #empty list of subs: ABCD,..etc
cyclength <- unlist(lapply(try.nice, nrow))
for(i in 1:10){ #running LSTAR model
  for(j in 2:5){
    #results for LSTAR model 
    lstarres[[i]][[j-1]] <- tryCatch({
      tsDyn::lstar(try.nice[[i]][,j], m=mdim, d=klag)}, #d = lag found through AIC
      error=function(e) list(fitted.values=rep(NA, (cyclength[[i]]-(klag*mdim))), 
                             residuals=rep(NA, (cyclength[[i]]-(klag*mdim))),
                             model.specific=list(coefficients = rep(NA, (4+2*mdim)))))
    #if error, output NA, (no fit) w/ fitted values and residual rep length times minus klag
  }
}
lstarres.fits <- lapply(lstarres, function(x) list(fits=x))
F = t(do.call(cbind, lapply(try.nice, function(x) x[-1][1:40,])))
meplz3 <- allresids(try.nice, modparams, lstarres.fits, Fer = F)

plot(x = 1:40, y = F[15,])


someone <- log(meplz)
