#general diagnostic 
model_comp <- function(subsets){
  
sig_est <- list(l4, b4, lstar)
res <- vector('list', length(sig_est))

for(sig in sig_est){
  
  if(sig %in% c(l4, b4)){
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
  if(identical(l5, sig)){
    res[[3]] <- resstat
      }
  if(identical(b5, sig)){
    res[[4]] <- resstat
      }
  }

if(sig == "lstar"){
    
  #running lstar function w/ certain lag
  lstarres <- vector("list", 8) #empty list of subs: ABCD,..etc
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
  #res1 is the list of subsets' fitted values and residuals
  res1 <- list()
  #matrix with fitted values and residuals as two columned vector
  for(i in 1:8){ #list of dataframes for each subset
    mat <- matrix( NA, nrow = 40, ncol = 25)
    mat[,1] <- 1:40 #different cycle lengths
    for(j in 1:12){
      vector <- c(lstarres[[i]][[j]]$fitted.values, lstarres[[i]][[j]]$residuals)
      ind1 = j*2 ; ind2 = ind1+1 #gets columns 2,3 ; 3,4; 5,6
      mat[(1+(klag*mdim)):40, ind1:ind2] <- matrix(vector, ncol=2) #fitted and residuals
    }
    res1[[i]] <- data.frame(mat) #df in dataframe
  }
  ###ADDED THRESHOLD CT VALUE: Beginning of RSSthres, RSS Red
  ###CT with Threshold STAR
  lstarcoef <- vector("list", sublength) #empty list of subs: ABCD,..etc
  cycCT.thover <- vector("list", sublength)
  cycCT.unl <- list() 
  #lstarcoef is list of lists of all coefficients 
  #cycCT.thover is list of lists of all cycles' LSTAR fitted > LSTAR threshold
  #cycCT.unl temporary unlisted CT.thover, because list of lists didn't work with which statement
  lstarres.unl = unlist(lstarres, recursive = FALSE)
  for(i in 1:10){
    for(j in 1:4){
      indicator <- j+4*(i-1)
      lstarcoef[[indicator]] <- lstarres.unl[[indicator]]$model.specific$coefficients #coeffs
      lstarcoef[[indicator]]["th"] <- abs(lstarcoef[[indicator]]["th"])
      if(sum(is.na(lstarcoef[[indicator]])) > 0){cycCT.unl[[indicator]] <- NA}
      else{
        #cycles greater than threshold in lagged format
        cycCT.unl[[indicator]] <- which(lstarres.unl[[indicator]]$fitted.values > lstarcoef[[indicator]]["th"])
      }
    }
    ind2 = 4*i ; ind1 = ind2-(4-1)
    cycCT.thover[[i]] <- cycCT.unl[ind1:ind2] #back to list of lists
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
  
    }
  }
}
