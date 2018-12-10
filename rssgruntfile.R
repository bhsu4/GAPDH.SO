#finding exponential phase with different methods
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/")
source("GAPDH.SO/genparams.R")
source("GAPDH.SO/rsslocal.R")

load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-520e_001119.Rda")
try.nice <- unlist.genparams(tst)
Fluo <- read_try(try.nice)

library(TTR) ; library(pracma)
exp_calc(method = 'RLE', subs = try.nice, thr=1.02) #RLE method
exp_calc(method = 'AQB', subs = try.nice) #AQB method
exp_calc(method = 'RAW', subs = try.nice) #RAW method
exp_calc(method = 'SSG', subs = try.nice) #SSG method
exp_calc(subs=try.nice, thr=1.02, all=TRUE) #all methods

for(i in 1:40){
  plot(Fluo[i,])
  abline(v=nicetestctau[i,2:3], lty='dashed', col='red')
  abline(v=nicetestctau[i,4:5], lty='dotted', col='green')
  abline(v=nicetestctau[i,6:7], lty='dotdash', col='blue')
  legend("topleft", c('rle', 'ema', 'mid'), 
         col=c('red', 'green', 'blue'), 
         lty=c('dashed', 'dotted', 'dotdash'))
}

#testing on another nice dataset
modparams <- lapply(try.nice, function(x) sub_genparams(listdf=x, est=l5)) #l5
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
meplz3 <- allresids(try.nice, modparams, lstarres.fits, Fer = Fluo)
