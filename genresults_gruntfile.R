#GAPDH.SO TESTING
#outside params generated
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

#GAPDH.SO results
gap.sig <- gen_results(subsets, modparams)
gap.lstar <- gen_results(subsets, lstarres.fits)

#miRcomp TESTING
#not nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst)
#nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams(tst)
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

#miRcomp results 
try.sig <- gen_results(try, modparams)
try.lstar <- gen_results(try, lstarres.fits)

#nice
modparams <- lapply(try.good, function(x) sub_genparams(listdf=x, est=sig))
lstarres <- vector("list", 10) #empty list of subs: ABCD,..etc
cyclength <- unlist(lapply(try.good, nrow))
for(i in 1:10){ #running LSTAR model
  for(j in 2:5){
    #results for LSTAR model 
    lstarres[[i]][[j-1]] <- tryCatch({
      tsDyn::lstar(try.good[[i]][,j], m=mdim, d=klag)}, #d = lag found through AIC
      error=function(e) list(fitted.values=rep(NA, (cyclength[[i]]-(klag*mdim))), 
                             residuals=rep(NA, (cyclength[[i]]-(klag*mdim))),
                             model.specific=list(coefficients = rep(NA, (4+2*mdim)))))
    #if error, output NA, (no fit) w/ fitted values and residual rep length times minus klag
  }
}
lstarres.fits <- lapply(lstarres, function(x) list(fits=x))

#miRcomp results 
trygood.sig <- gen_results(try.good, modparams)
trygood.lstar <- gen_results(try.good, lstarres.fits)
