gen_results <- function(subs, params){
  ##define lengths
  sublength = length(subs)
  repnames <- unique(unlist(lapply(subsets, names)))[!grepl("Cycle", unique(unlist(lapply(subsets, names))))]
  replength = length(repnames)
  ##finding resids 
  modresids <- lapply(params, function(x) lapply(x$fits, resid))
  ##unlisting resids
  modresids.unl <- unlist(modresids, recursive = FALSE)
  #lengths of the resids, cycle lengths
  modresids.length <- lapply(modresids.unl, length)
  #finding rss total
  modrss <- lapply(modresids, function(x) lapply(x, function(x) sum(x[which(!is.na(x))]^2)))
  
  ##autocorrelation and dw for amplification and residual values
  reg.amp <- vector('list', sublength) ; reg.res <- vector('list', sublength) #dynlm formula results
  dw.amp <- list() ; dw.res <- vector('list', sublength) #dw-test values
  for(i in 1:sublength){
    reg.amp[[i]] <- lapply(subs[[i]][ ,2:length(subs[[i]])], 
                           function(y) dynlm(y ~ subs[[i]][,1])) #all amp can run
    reg.res[[i]] <- lapply(modresids[[i]], function(y) try(dynlm(y ~ subs[[i]][,1]), silent=TRUE))
    #residuals contain NA
    dw.amp[[i]] <- lapply(reg.amp[[i]], function(x) durbinWatsonTest(x)) #amplification values
  }
  for(i in 1:sublength){ #residuals
    for(j in 1:(length(subs[[i]])-1)){ #NA resids cannot run
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
  dw.res <- lapply(dw.res, function(x) lapply(x, function(x) x[1:3]))
  dw.resmat <- matrix(unlist(dw.res), ncol = 3, byrow=TRUE)
  #one dw
  dw.res <- cbind(dw.ampmat, dw.resmat)
  
  ##CT values: sigmoidal and lstar  
  if(isTRUE(params[[1]]$fits[[1]]$MODEL$name %in% c("l5", "l4", "b5", "b4")) == TRUE){
    #finding CT value
    mod1 <- lapply(subsets, function(x) modlist(x, model = l5))
    res1 <- lapply(mod1, function(x) getPar(x, type = "curve", cp = "cpD2", eff = "sliwin")) 
    res1CT <- lapply(res1, function(x) x[1,]) #first row is CT
    #specific CT +/- 2 cycles
    rssgrey <-list()
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        indij = ((replength/sublength)*(i-1))+j  
        uppercyc = floor(unlist(res1CT)[[indij]]+2)  #+2 cycles
        #klag-adjusted lower cycle
        lowercyc = ceiling(unlist(res1CT)[[indij]]-2) #-2cycles }
        #list of lists of rss grey region
        rssgrey[[indij]] <- sum(modresids.unl[[indij]][lowercyc:uppercyc]^2) 
      }
    }
    #rss.mat for sigmoidal
    rss.mat <- matrix(c(matrix(unlist(modrss), ncol=1), 
                        matrix(unlist(rssgrey), ncol=1), 
                        matrix(unlist(res1CT), ncol=1)), ncol=3)
  }
  #LSTAR model
  else if(is.null(params[[1]]$fits[[1]]$MODEL$name) == TRUE){
    ###CT with Threshold STAR
    cycCT.thover <- vector("list", sublength)
    #lstarcoef is list of lists of all coefficients 
    #cycCT.thover is list of lists of all cycles' LSTAR fitted > LSTAR threshold
    #cycCT.unl temporary unlisted CT.thover, because list of lists didn't work with which statement
    lstarres.unl = unlist(lstarres.fits, recursive = FALSE)
    lstarcoef <- lapply(lstarres.unl, function(x) lapply(x, function(y) y$model.specific$coefficients))
    lstarcoef.th <- lapply(lstarcoef, function(x) lapply(x, function(y) y["th"]))
    
    #navalues <- lapply(lstarcoef.th, is.na)
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        if(sum(is.na(lstarcoef[[i]][[j]])) > 0){cycCT.unl[[i]][[j]] <- NA}
        else{
          #cycles greater than threshold in lagged format
          cycCT.thover[[i]][[j]] <- which(lstarres.unl[[i]][[j]]$fitted.values > lstarcoef.th[[i]][[j]])
        }
      }
    }
    overth.diff <- function(x){ #diff takes lagged difference, rle shows lengths of lagged diff
      #ex. x = 4,5,6,7,10,13,25,26,27,28,29,30
      #lengths are 3,2,1,5 these are how many of those differences (max this)
      #values are 1,3,12,1 these are diffrences (this needs to be = 1)
      if(sum(is.na(x)) == 0){ rle(diff(x)) }
      else{ list(lengths=NA) }
    }      
    overth.diffcyc <- function(x){ #fitting conditions above
      #max lengths and choose values of 1
      #logic here is that: in cases where more fluctuations, more noise, multiple times when
      #cycles cross threshold. take the cycle that has most cycles after threshold that is 
      #higher than the threshold continuously
      if(sum(is.na(x)) == 0){ which(x$lengths == max(x$lengths) & x$values == 1) }
      else{ NA }
    }
    #difftimes are rle (run length encoding) lagged diff
    difftimes <- lapply(cycCT.thover, function(x) lapply(x, overth.diff))
    #cycle index for cycle in list of those > threshold
    index.length <- lapply(difftimes, function(x) lapply(x, overth.diffcyc))
    index.unl <- unlist(index.length, recursive = FALSE)
    index.unl <- lapply(index.unl, max) #if two same length of lagged diff, choose latter cycles
    #final CT list
    cycCT.th <- vector("list", sublength)
    #find the intermediate CT value (intermediate means 0.5) surpassing threshold
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        ind = j+(replength/sublength)*(i-1)
        if(sum(is.na(index.unl[[ind]])) > 0){ #any NAs get NA
          cycCT.th[[i]][[j]] <- NA
        }
        else if(index.unl[[ind]] > 1){ #more than one string of consec > threshold
          #ex. 1,2,3,4,5,6,7,8,9,10,.....,24 (20th term), 28,29,30,....
          #difftimes lengths: 11, 1, 5, 1, 1, 1, 11
          #          values:  1, 3, 1, 3, 1, 4, 1
          #index.unl = 7, take sum of 1:6 of the lengths = 20, so 21st is our chosen pt
          tmd <- sum(difftimes[[i]][[j]]$lengths[1:(index.unl[[ind]]-1)]) #sum of before
          cycCT.th[[i]][[j]] <- (cycCT.thover[[i]][[j]][[tmd+1]] + (cycCT.thover[[i]][[j]][[tmd+1]]-1))/2 + (mdim*klag) #actual pt after add
        }
        else{ 
          #at index = 1, just take that first point, and point-1 , find average
          cycCT.th[[i]][[j]] <- (cycCT.thover[[i]][[j]][[1]] + (cycCT.thover[[i]][[j]][[1]]-1))/2 + (mdim*klag) #actual pt after addition of mdim*klag
        }
      }
    }
    ##rsslstarr THRESHOLD CT for +/- 2 cycles
    #cycles for CT LSTAR THRESHOLD
    unlistcycCTthr <- unlist(unlist(cycCT.th, recursive=FALSE))
    #add NA to those that are outside lower bound first
    unlistcycCTthr[which(unlistcycCTthr < (0.5+klag*mdim+2))] <- NA
    #rss for red box
    rsslstarr <- vector('list', sublength) 
    for(i in 1:(sublength)){
      for(j in 1:(replength/sublength)){
        #combined indicator with i and j
        indij = ((replength/sublength)*(i-1))+j  
        #all NA values
        if(sum(is.na(unlistcycCTthr[[indij]])) > 0){
          #this does not go through, will have to run this at the end
          rsslstarr[[i]][[j]] <- NA
        }
        #filter out low values, change to NA unlistcycCTthr
        else if(sum(isTRUE(unlistcycCTthr[[indij]] > (cyclength[[i]]-(klag*mdim+1)))) > 0){ #strict upper bound
          uppercyc = cyclength[[i]]-(klag*mdim)
          lowercyc = ceiling(unlistcycCTthr[[indij]]-(klag*mdim)-2)}
        else{    #(sum(is.na(unlistcycCTthr[[indij]])) == 0){ 
          uppercyc = floor(unlistcycCTthr[[indij]]-(klag*mdim)+2)  #+2 cycles
          #klag-adjusted lower cycle
          lowercyc = ceiling(unlistcycCTthr[[indij]]-(klag*mdim)-2)} #-2cycles 
        #list of lists of rss red region
        rsslstarr[[i]][[j]] <- sum(lstarres[[i]][[j]]$residuals[lowercyc:uppercyc]^2)
        #rsslstarr provies RSS for NA values, we undo this here (rewrite above, need both)
        if(sum(is.na(unlistcycCTthr[[indij]])) > 0){ 
          rsslstarr[[i]][[j]] <- NA 
          print(paste(k, "-", indij, "no CTth"))}
        else{}
      }
    }
    #unlisted 1:(replength) of rss grey region
    unlist.rsslstarr= unlist(rsslstarr, recursive = FALSE)
    #rss.mat matrix for LSTAR
    rss.mat <- matrix(c(matrix(unlist(modrss), ncol=1), unlist.rsslstarr, unlistcycCTthr), ncol=3)
  } #else if lstar
  else{ print("no model statistics") }
  
  ##combining all into matrix
  empmat <- matrix(data=NA, nrow=replength, ncol=10)
  empmat[,1] <- unique(unlist(lapply(subsets, colnames)))[-grep("Cycle", unique(unlist(lapply(subsets, colnames))))]
  empmat[,2:4] <- matrix(rss.mat, ncol=3)
  empmat[,5:10] <- dw.res
  #dataframe output
  resstat <- data.frame(empmat)
  colnames(resstat) <- c("Replicate", "RSS", "RSSlocal", "CT", 
                         "amp.r", "amp.dw", "amp.p", 
                         "res.r", "res.dw", "res.p")   
  return(resstat)
}

#outside params generated
modparams <- lapply(subsets, function(x) sub_genparams(listdf=x, est=sig))
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

#GAPDH.SO
gap.sig <- gen_results(subsets, modparams)
gap.lstar <- gen_results(subsets, lstarres.fits)
