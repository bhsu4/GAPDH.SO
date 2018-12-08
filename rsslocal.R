branchfunc <- function(subs, method='RLE', plot=FALSE){
  #setting lengths
  sublength = length(subs)
  repnames <- unique(unlist(lapply(subs, names)))[!grepl("Cycle", unique(unlist(lapply(subs, names))))]
  replength = length(repnames)
  cyclength <- unlist(lapply(subs, nrow))
  #creating Fluo 
  max_len <- max(unlist(lapply(subs, function(x) lapply(x, length))))
  subs.unl <- unlist(lapply(subs, function(x) x[-1]), recursive = FALSE)
  empmat <- matrix(NA, max_len, replength) #set max cycles
  for(i in 1:40){ empmat[1:length(subs.unl[[i]]), i] <- subs.unl[[i]] } #NA filled for non-max lengths
  Fluo = t(data.frame(empmat)) #transposed 
  rownames(Fluo) <- unlist(lapply(LETTERS[1:10], function(x) paste0(x, 1:4))) #names assigned
  
  #call exponential calculation
  exp_ctau <- exp_calc(method, subs)
  ctau1 <- unlist(exp_ctau[,2])
  ctau2 <- unlist(exp_ctau[,3])
  
  #finding initial
  yn <- matrix(0, r, 1)
  yn_1 <- matrix(0, r, 1)
  p_tilde <- matrix(0, r, 1)
  mu_b <- matrix(0, r, 1)
  var_mu_b <- matrix(0, r, 1)
  for (i in 1:r) {
    yn[i] <- sum(Fluo[i, ctau1[i]:ctau2[i]], na.rm=TRUE)
    yn_1[i] <- sum(Fluo[i, ctau1[i]:(ctau2[i] - 1)])
    p_tilde[i] <- (yn[i] - Fluo[i, ctau1[i]] - yn_1[i])/yn_1[i]
    mu_b[i] <- p_tilde[i]/((1 + p_tilde[i])^n) * yn[i] * (1 - (1 + p_tilde[i])^(ctau1[i] - n))^(-1)
  }
  n_an <- matrix(0, r, n) ; cdivas <- matrix(0, r, 1) ; appF <- matrix(0, r, n)
  n_an[,1] <- mean(mu_b)
  for(i in 1:r){
    for(j in 2:n){
      n_an[i,j] <- n_an[i, j-1] + rbinom(1, round(n_an[i, j-1]), p_tilde[i])
    }
    cdivas[i,] = (n_an[i, round((ctau1[i]+ctau2[i])/2)]/Fluo[i, round((ctau1[i]+ctau2[i])/2)])
    appF[i,] = n_an[i,]/cdivas[i,]
  }
  #backtransform constant
  backconst <- matrix(0, r, n) ; res_aq <- matrix(0, r, n)
  for(i in 1:r){
    backconst[i,] = cumsum(Fluo[i,])/Fluo[i,]
    res_aq[i,] <- cumsum(appF[i,])/backconst[i,]
  }
  #plotting
  if(plot){  
    for(i in 1:r){
      plot(x=1:max_len, y=res_aq[i,], col=2)
      points(x=1:max_len, y=Fluo[i,])
    }
  }
  #resid.res list of list
  resid.res <- vector('list', sublength)
  for(i in 1:sublength){
    for(j in 1:(replength/sublength)){
      indij = (replength/sublength)*(i-1)+j
      resid.res[[i]][[j]] <- res_aq[indij,] - Fluo[indij,]
    }
  }
  return(list(fitted.values = res_aq, residuals = resid.res, expo = exp_ctau))
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
  branch.res <- branchfunc(subs)
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
try.good <- unlist.genparams(tst)
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
meplz3 <- allresids(try.good, modparams.good, lstarres.fits, Fer = Fluo)


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

sdres <- apply(meplz, 2, mean, na.rm = TRUE)
sdres2 <- apply(meplz, 2, sd, na.rm = TRUE)


###finding the exponential
read_try <- function(subs){
  #compute lengths first
  max_len <- max(unlist(lapply(subs, function(x) lapply(x, length))))
  subs.unl <- unlist(lapply(subs, function(x) x[-1]), recursive = FALSE)
  reps_len <- length(subs.unl) #40
  subs_len <- length(subs) #10
  empmat <- matrix(NA, max_len, reps_len)
  for(i in 1:reps_len){
    empmat[1:length(subs.unl[[i]]), i] <- subs.unl[[i]]
  }
  Fluo = t(data.frame(empmat)) #row-wise fluorescence
  a1 <- dim(Fluo) ; r <- a1[1] ; n <- a1[2] #row-wise fluo dimensions
  rownames(Fluo) <- unlist(lapply(LETTERS[1:subs_len], function(x) paste0(x, 1:(reps_len/subs_len))))
  return(Fluo)
}

exp_calc <- function(method = 'RLE', subs, thr=1.02, all = FALSE){
  
  #compute lengths first
  max_len <- max(unlist(lapply(subs, function(x) lapply(x, length))))
  subs.unl <- unlist(lapply(subs, function(x) x[-1]), recursive = FALSE)
  reps_len <- length(subs.unl) #40
  subs_len <- length(subs) #10
  empmat <- matrix(NA, max_len, reps_len)
#  for(i in 1:reps_len){
#    empmat[1:length(subs.unl[[i]]), i] <- subs.unl[[i]]
#  }
#  Fluo = t(data.frame(empmat)) #row-wise fluorescence
  a1 <- dim(Fluo) ; r <- a1[1] ; n <- a1[2] #row-wise fluo dimensions
#  rownames(Fluo) <- unlist(lapply(LETTERS[1:subs_len], function(x) paste0(x, 1:(reps_len/subs_len))))
  len_used <- matrix(unlist(lapply(subs.unl, length)), ncol=1)
  Fluo <- read_try(subs)
    
  #run-length encoding method
  if (method == 'RLE' & all == FALSE || all == TRUE){
    ctau1_rle <- matrix(NA, reps_len, 1) ; ctau2_rle <- matrix(NA, reps_len, 1)
    for(i in 1:reps_len){
      eff_est <- rle( (Fluo[i, 2:n]/Fluo[i, 1:(n - 1)]) > thr ) #run length encoding
      tester.mat <- matrix(eff_est$lengths, nrow=1) 
      overth <- which(eff_est$values == TRUE) #all indexes over thr
      max.overth <- max(eff_est$lengths[overth]) #max consec cycles at the indexes
      ctau1_rle[i,] <- sum(tester.mat[,1:(which(eff_est$lengths == max.overth & eff_est$values == TRUE)-1)])+1
      ctau2_rle[i,] <- ctau1_rle[i,] + max.overth
    }
    
    #ema with RLE
    ctau2_ema <- matrix(NA, reps_len, 1)
    for(i in 1:reps_len){
      Fluo_EMA3 <- EMA( (Fluo[i, 2:n]/Fluo[i, 1:(n-1)])[ctau1_rle[i]:ctau2_rle[i]-1] , n=3) #window=3
      thr_ema <- mean(Fluo_EMA3[which(Fluo_EMA3>0)]) #mean of Fluo-EMA3 > 0
      ctau2_ema[i,] <- ctau1_rle[i] + max(which(Fluo_EMA3 > thr_ema))
    }
    
    #midpoint method
    ctau1_mid <- matrix(NA, reps_len, 1) ; ctau2_mid <- matrix(NA, reps_len, 1)
    #this works, consider chaning threshold 1.25 to something else
    for(i in 1:reps_len){
      midpt <- round((ctau1_rle[i]+ctau2_rle[i])/2)
      range1 <- midpt-1 ; range2 <- midpt+1 #starting 3-cycle range
      Feff <- Fluo[i, 2:n]/Fluo[i, 1:(n-1)]
      while(range1 > ctau1_rle[i] & range2 < ctau2_rle[i]){
        test1 <- range1-1 ; test2 <- range2+1 
        comp_mean1 <- mean(Feff[test1:range2])
        comp_mean2 <- mean(Feff[range1:test2])
        #print(paste(comp_mean1, comp_mean2))
        #finding larger mean range
        if(comp_mean1 > comp_mean2){
          range1 = test1 ; range2 = range2 
          if(comp_mean1 <= 1.25){ break }
        }
        else if(comp_mean1 < comp_mean2){ 
          range1 = range1 
          range2 = test2 
          if(comp_mean2 <= 1.25){ break }
        }
        #print(paste(i, mean(Feff[range1:range2]), range1, range2))
      }
      #print(paste(i, range1, range2))
      #cat("\n")
      ctau1_mid[i,] <- range1 ; ctau2_mid[i,] <- range2
    }
    
    ctau_rle <- data.frame(cbind(len_used, ctau1_rle, ctau2_rle, 
                                 ctau1_rle, ctau2_ema, 
                                 ctau1_mid, ctau2_mid))
    colnames(ctau_rle) <- c("len", "rle1", "rle2", "ema1", "ema2", "mid1", "mid2")
    if(all == FALSE){ return(ctau_rle) }
  }
  
  #AQB method
  else if (method == 'AQB' & all == FALSE || all == TRUE){
    tau1 <- 1.02 ; tau2 <- 1.02
    ctau1_aqb <- matrix(0, r, 1) ; ctau2_aqb <- matrix(0, r, 1)
    for (i in 1:r) {
      ctau1_aqb[i] <- min(which(Fluo[i, 2:n]/Fluo[i, 1:(n - 1)] >= tau1))
      if (any(Fluo[i, (ctau1_aqb[i] + 1):n]/Fluo[i, ctau1_aqb[i]:(n - 1)] <= tau2)) {
        ctau2_aqb[i] <- min(which(Fluo[i, (ctau1_aqb[i] + 1):n]/Fluo[i, ctau1_aqb[i]:(n - 1)] <= tau2)) + ctau1_aqb[i]
      }
      else {
        ctau2_aqb[i] = n
      }
    }
    ctau_aqb <- data.frame(cbind(len_used, ctau1_aqb, ctau2_aqb))
    colnames(ctau_aqb) <- c("len", "aqb1", "aqb2")
    if(all == FALSE){ return(ctau_aqb) }
  }
  
  #RAW method: using efficieincy and slope of raw fluorescence
  else if (method == 'RAW' & all == FALSE || all == TRUE){
    #cumulative averages
    cumulative_left <- matrix(NA, reps_len, max_len) ; threshold <- matrix(NA, reps_len, 1)
    ctau1_raw <- matrix(NA, reps_len, 1) ; ctau2_raw <- matrix(NA, reps_len, 1)
    for(i in 1:reps_len){
      subs_len <- length(Fluo[i,][which(!is.na(Fluo[i,]))]) #cycle length
      threshold[i,] = 0.04 * (max(Fluo[i,], na.rm=TRUE) - min(Fluo[i,], na.rm=TRUE))
      for(n in 1:max_len){
        cumulative_left[i,n] = sum(Fluo[i,1:n])/sum(!is.na(Fluo[i,1:n]))
      }
      ctau1_raw[i,] <- min(which(Fluo[i,] - cumulative_left[i,] > threshold[i,]))
      ctau2_raw[i,] <- ctau1_raw[i,] + which.max((Fluo[i, 2:n] - Fluo[i, 1:(n-1)])[ctau1_raw[1,]:(subs_len)])+1
    }
    ctau_raw <- data.frame(cbind(len_used, ctau1_raw, ctau2_raw))
    colnames(ctau_raw) <- c("len", "raw1", "raw2")
    if(all == FALSE){ return(ctau_raw) }
  }
  
  #SSG method: savitzsky-golay smoothing
  else if (method == 'SSG' & all == FALSE || all == TRUE){
    #savitzky-golay smoothing local min/max method
    ctau1_ssg <- matrix(NA, reps_len, 1) ; ctau2_ssg <- matrix(NA, reps_len, 1)
    for(i in 1:reps_len){
      subs = Fluo[i,][which(!is.na(Fluo[i,]))]
      subs_len <- length(subs)
      if(sd(subs) < 100){ #noisy data
        point_filt = savgol(subs, fl=25, forder=3) #window=25, poly-order=3
        sslope_filt = savgol((point_filt[2:(subs_len)]-point_filt[1:(subs_len-1)]), fl=15, forder=3)
      } else{ #clean data
        point_filt = savgol(subs, fl=3, forder=3) #window=5, poly-order=3
        sslope_filt = savgol((point_filt[2:(subs_len)]-point_filt[1:(subs_len-1)]), fl=5, forder=3)
      }
      #restandardize
      rs.point_filt = (point_filt - min(point_filt))/(max(point_filt) - min(point_filt))
      rs.sslope_filt = (sslope_filt - min(sslope_filt))/(max(sslope_filt) - min(sslope_filt))
      #find cycles
      lmin_val = localMinima(rs.sslope_filt) #lmin cyc
      lmax_val = localMaxima(rs.sslope_filt) #lmax cyc
      cyc_lmax = which(rs.sslope_filt == 1) #standardized max at 1
      cyc_lmin = lmin_val[max(which(lmin_val < cyc_lmax))] #first lmin value < lmax
      #output
      ctau1_ssg[i,] <- cyc_lmin + 1
      ctau2_ssg[i,] <- cyc_lmax + 1
    }
    ctau_ssg <- data.frame(cbind(len_used, ctau1_ssg, ctau2_ssg))
    colnames(ctau_ssg) <- c("len", "ssg1", "ssg2")
    if(all == FALSE){ return(ctau_ssg) }
  }
  #final return output
  if(all == TRUE){ 
    ctau_all = data.frame(cbind(len_used, ctau1_rle, ctau2_rle,  #rle
                                ctau1_rle, ctau2_ema,  #ema
                                ctau1_mid, ctau2_mid,  #mid-seq
                                ctau1_aqb, ctau2_aqb,  #aqb
                                ctau1_raw, ctau2_raw,  #raw
                                ctau1_ssg, ctau2_ssg)) #savit-golay
    colnames(ctau_all) <- c("len", "rle1", "rle2", "ema1", "ema2", "mid1", "mid2",
                            "aqb1", "aqb2", "raw1", "raw2", "ssg1", "ssg2")
    return(ctau_all)
  }
}

localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}




