#finding exponential phase with different methods

load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-520e_001119.Rda")
try.nice <- unlist.genparams(tst)
Fluo <- read_try(try.nice)

library(TTR) ; library(pracma)
exp_calc(method = 'RLE', subs = try.good, thr=1.02) #RLE method
exp_calc(method = 'AQB', subs = try.good) #AQB method
exp_calc(method = 'RAW', subs = try.good) #RAW method
exp_calc(method = 'SSG', subs = try.good) #SSG method
exp_calc(subs=try.good, thr=1.02, all=TRUE) #all methods

for(i in 1:40){
  plot(Fluo[i,])
  abline(v=nicetestctau[i,2:3], lty='dashed', col='red')
  abline(v=nicetestctau[i,4:5], lty='dotted', col='green')
  abline(v=nicetestctau[i,6:7], lty='dotdash', col='blue')
  legend("topleft", c('rle', 'ema', 'mid'), 
         col=c('red', 'green', 'blue'), 
         lty=c('dashed', 'dotted', 'dotdash'))
}