library(miRcompData)
data("miRcompData")
miRcompData2 <- miRcompData
miRcompData2 <- miRcompData2[!(miRcompData2$SampleID == "E1_" | miRcompData2$SampleID == "E2_"),] #delete E's

#*difference between subset and whole is that subset using only 1 fluoresence and 
#*whole uses a representative group fluorescence

#using targlist
source("GAPDH.SO/targlist.R")
#singtarget.list for pulling up 1 target
targetatt <- singtarget.list(miRcompData2, target = targnames[1]) 
#saving all targets
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
savetarget.list(miRcompData2) #saves all the miRcompData as list of lists

#generating parameter estimates and residuals (whole + subset)
source("GAPDH.SO/genparams.R")
#after all targets -- loads as tst, replaces
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst) #list of lists organized here
KW1_tstresults <- sub_genparams(l4, try[[1]]) #specfic subset results
subplot_resid(try[[1]], KW1_tstresults) #specific subset residuals
result.try <- genparams(b5, try) #whole target params
plot_resid(try, result.try) #whole target resids

#plotting amplification surve + residual plot
##DW-statistic package loading
library(dynlm) ; library(car)

#not nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst)
#nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams(tst)

#plotting fluo curves + finding parameter est in matrix
source("GAPDH.SO/plot_sig.R")
plot_sig(l4, try, plot=TRUE) #plots shown
plot_sig(b5, try) #no plots shown
plot_sig(l5, try, macro=1, z=4, plot=TRUE) #tester for NA point (pulls out try[[4]][[1]])
plot_sig(b5, try, macro=3, z=4, plot=TRUE) #tester for good point
plot_sig(l4, try, z=4, plot=TRUE)

#smaller data set with only 4 targets
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello")
library(foreach) ; library(dynlm) ; library(car) ; library(profvis) #foreach for parallel execution
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello", 
                pattern =  "^targ_")
profvis({sig_est(l5, miRcompData2, getfiles=getfiles) }) #system periodicity time estimation

#larger data set
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont", 
                pattern =  "^targ_")

#running to find matrix combined
l5dat <- sig_est(l5, miRcompData2, getfiles=getfiles)
l4dat <- sig_est(l4, miRcompData2, getfiles=getfiles)
b5dat <- sig_est(b5, miRcompData2, getfiles=getfiles)
b4dat <- sig_est(b4, miRcompData2, getfiles=getfiles)

#saving files
save(l5dat, file = paste0("l5_" , "dat", ".Rda"))
save(l4dat, file = paste0("l4_" , "dat", ".Rda"))
save(b5dat, file = paste0("b5_" , "dat", ".Rda"))
save(b4dat, file = paste0("b4_" , "dat", ".Rda"))

#loading in matrices of parameter est
getmat <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont", 
              pattern = "_dat")

load(file = getmat[[1]]) ; load(file = getmat[[2]]) ; load(file = getmat[[3]]) ; load(file = getmat[[4]])
l5dat = help1 ; l4dat = help2 ; b5dat = help3 ; b4dat = help4

