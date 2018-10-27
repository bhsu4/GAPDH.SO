#outside params generated
modparams <- lapply(subsets, function(x) sub_genparams(listdf=x, est=sig))
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

lstarres.fits <- lapply(lstarres, function(x) list(fits=x))
#finding resids same
modresids <- lapply(modparams, function(x) lapply(x$fits, resid))
lstarresresids <- lapply(lstarres.fits, function(x) lapply(x$fits, resid))
#unlisting same
modresids.unl <- unlist(modresids, recursive = FALSE)
modresids.length <- lapply(modresids.unl, length)
lstarresresids.unl <- unlist(lstarresresids, recursive = FALSE)
lstarresresids.length <- lapply(lstarresresids.unl, length)
#finding rss total
rss_na <- function(x){
  x <- x[which(!is.na(x))]
  sum(x^2)
}
modrss <- lapply(modresids, function(x) lapply(x, function(x) sum(x[which(!is.na(x))]^2)))
lstarrss <- lapply(lstarresresids, function(x) lapply(x, function(x) sum(x[which(!is.na(x))]^2)))

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

#####for lstar model same for all 
reg.amp <- vector('list', 8) ; reg.res <- vector('list', 8) #dynlm formula results
dw.amp <- list() ; dw.res <- vector('list', 8) #dw-test values
for(i in 1:8){
  reg.amp[[i]] <- lapply(subsets[[i]][ ,2:length(subsets[[i]])], 
                         function(y) dynlm(y ~ subsets[[i]][,1])) #all amp can run
  reg.res[[i]] <- lapply(lstarresresids[[i]], function(y) try(dynlm(y ~ subsets[[i]][,1]), silent=TRUE))
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



###CT value extraction is different 
modparams$
