library(miRcompData)
data("miRcompData")
miRcompData2 <- miRcompData
miRcompData2 <- miRcompData2[!(miRcompData2$SampleID == "E1_" | miRcompData2$SampleID == "E2_"),]
targnames <- unique(c(miRcompData$TargetName)) 
sampnames <- unique(c(miRcompData$SampleID))

source("GAPDH.SO/targlist.R")
targetatt <- singtarget.list(miRcompData2, target = targnames[1])
savetarget.list(miRcompData2) #saves all the miRcompData as list of lists

setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
savetarget.list(miRcompData2) #saves all the miRcompData as list of lists

source("GAPDH.SO/genparams.R")
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst)
result.try <- genparams(b5, try)
plot_resid(try, result.try)
KW1_tstresults <- sub_genparams(l4, try[[1]])
subplot_resid(try[[1]], KW1_tstresults)

#running through all files (??????)
files <- list.files(path="C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets", 
                    pattern="*.Rda", full.names=T, recursive=FALSE)
lapply(files, function(x) {
  t <- read.table(x, header=T) # load file
  # apply function
  out <- function(t)
    # write to file
    write.table(out, "path/to/output", sep="\t", quote=F, row.names=F, col.names=T)
})

##DW-statistic
library(dynlm)
library(car)
reg1 = dynlm(A1 ~ Cycle, data = try[[1]])
hello <- durbinWatsonTest(reg1)

tst <- matrix( ,nrow = 10, ncol = 4)
for(i in 1:10){
  for(k in 1:4){
    reg <- dynlm(try[[i]][[k]] ~ Cycle, data = try[[i]])
    hello <- durbinWatsonTest(reg)qq
    tst[i,k] <- hello$p
  }
}
plot(x=1:10, y=tst[,1], xlim=c(1,10), col=1, ylim=range(tst))
for(i in 2:4){
  points(x=1:10, y=tst[,i], col = i)
}

##CT Value
ml1 <- modlist(try[[1]], model = l4)
res1 <- getPar(ml1, type = "curve", cp = "cpD2", eff = "sliwin")
str(res1)
res1[,1:4]

#mccall commands
files <- fir("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/", 
             pattern = "Rda$")
fff <- function(fname){
  load(paste0("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets", fname))
  print(str(tst))
}
fff(files[10])
sapply(files[1:10], fff)
fff(files[10])

##fitting the curves
source("GAPDH.SO/b5_model.R")
plot_b5(try, result.try)

#fix b5 to l4 / change every model
KW1_tstresults2 <- sub_genparams(l4, try[[1]])
subplot_l4(try[[1]], KW1_tstresults2)

KW_tstresults <- lapply(try, function(x) sub_genparams(l4, x))
for(i in 1:10){
  subplot_l4(try[[i]], KW_tstresults[[i]])
}

#nosignal 
par(mfrow=c(2,1))
par(oma=c(4,4,4,4),mar=c(0.25,0.25,0,0))

subsubplot_l4 <- function(listdf, par, k){
  xs = listdf$Cycle
  plot(x=xs, y=l4_model(xs, b=par$params$b[k], c=par$params$c[k],
                            d=par$params$d[k], e=par$params$e[k]), type="l",  
       xlab="Cycle", ylab="Fluorescence", 
       ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                    names(unlist(listdf)))])])), col = k, 
       xaxt = "n")
  points(x=xs, y=listdf[,k+1], cex=0.45)
  legend("topleft", c(names(listdf)[k+1]), 
         col=k, lty=1, cex=0.65)
  #mtext(paste("b=", round(par$params$b[k], 3)), side=3, adj=1, line=3)
}
subsubplot_l4(try[[1]], KW1_tstresults2, 3)

subsubplot_resid <- function(listdf, params, k){
    resids <- lapply(params$fits, resid)
    plot(y = resids[[k]][1:length(listdf$Cycle)], 
         x = params$fits[[k]]$DATA$Cycles[1:length(listdf$Cycle)], 
         ylim = range(unlist(resids)), 
         xlab = "Cycle", ylab="Fluorescence")
    abline(h=0)
    rss <- sum(resids[[k]]^2)
    
    reg <- dynlm(resids[[k]] ~ listdf$Cycle)
    hello <- durbinWatsonTest(reg)
    tog <- list(hello, rss)
    return(tog)
}
subsubplot_resid(try[[1]], KW1_tstresults, 3)

result.try <- genparams(l4, try)

reg <- dynlm(resids.tst[[3]] ~ try[[1]]$Cycle)

values <- data.frame(apply(result.try$params[3,], c(1,2), as.numeric))
values[,5:7] <- c(AC = hello$r, DW = hello$dw, p = hello$p)

subsubplots_tog <- function(listdf, par, k, res){
  par(mfrow=c(2,1))
  par(oma=c(4,4,4,4),mar=c(0.25,0.25,0,0))
  
  subsubplot_l4(listdf, par, k)
  tog <- subsubplot_resid(listdf, par, k) 
  values <- data.frame(apply(res$params[k,], c(1,2), as.numeric))
  values[,(length(values)+1):(length(values)+4)] <- c(tog[[1]]$r, 
                                                      tog[[1]]$dw, 
                                                      tog[[1]]$p,
                                                      tog[[2]])
  names(values) <- c(names(res$params), c("r", "dw", "p", "rss"))
  return(values)
}
nosig1_KW1_A3 <- subsubplots_tog(try[[1]], KW1_tstresults, 3, result.try)


##decent fits
KW8_tstresults <- sub_genparams(l4, try[[8]])
subplot_l4(try[[8]], KW8_tstresults)
subplot_resid(try[[8]], KW8_tstresults)

sig_KW8_H4 <- subsubplots_tog(try[[8]], KW8_tstresults, 4, KW10_tstresults)

KW10_tstresults <- sub_genparams(l4, try[[10]])
subplot_l4(try[[10]], KW10_tstresults)
subplot_resid(try[[10]], KW10_tstresults)

sig_KW10_J3 <- subsubplots_tog(try[[10]], KW10_tstresults, 3, KW10_tstresults)

#poor fits
KW5_tstresults <- sub_genparams(l4, try[[5]])
subplot_l4(try[[5]], KW5_tstresults)
subplot_resid(try[[5]], KW5_tstresults)

sig_KW5_E4 <- subsubplots_tog(try[[5]], KW5_tstresults, 4, KW5_tstresults)












library(dynlm) ; library(car)
source("GAPDH.SO/genparams.R")

###everything as function
subsubplots_tog <- function(est, listdf, k){
  #genearting parameters of specific
  par <- sub_genparams(est, listdf)
  #two graphs on top each other
  par(mfrow=c(2,1))
  par(oma=c(4,4,4,4),mar=c(0.25,0.25,0,0))
  #plotting amplification curve
  xs = listdf$Cycle
  if(est$name == "l4"){
  plot(x=xs, y=l4_model(xs, b=par$params$b[k], c=par$params$c[k],
                            d=par$params$d[k], e=par$params$e[k]), type="l",  
       xlab="Cycle", ylab="Fluorescence", 
       ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                    names(unlist(listdf)))])])), col = k, xaxt = "n")
  points(x=xs, y=listdf[,k+1], cex=0.45)
  legend("topleft", c(names(listdf)[k+1]), 
         col=k, lty=1, cex=0.65)
  }
  if(est$name == "b4"){
  plot(x=xs, y=b4_model(xs, b=par$params$b[k], c=par$params$c[k],
                            d=par$params$d[k], e=par$params$e[k]), type="l",  
       xlab="Cycle", ylab="Fluorescence", 
       ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                    names(unlist(listdf)))])])), col = k, xaxt = "n")
  points(x=xs, y=listdf[,k+1], cex=0.45)
  legend("topleft", c(names(listdf)[k+1]), 
         col=k, lty=1, cex=0.65)
  }
  if(est$name == "l5"){
    plot(x=xs, y=l5_model(xs, b=par$params$b[k], c=par$params$c[k],
                              d=par$params$d[k], e=par$params$e[k],
                              f=par$params$f[k]), type="l",  
         xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                      names(unlist(listdf)))])])), col = k, xaxt = "n")
    points(x=xs, y=listdf[,k+1], cex=0.45)
    legend("topleft", c(names(listdf)[k+1]), 
           col=k, lty=1, cex=0.65)
  }
  if(est$name == "b5"){ 
    plot(x=xs, y=b5_model(xs, b=par$params$b[k], c=par$params$c[k],
                              d=par$params$d[k], e=par$params$e[k],
                              f=par$params$f[k]), type="l",  
         xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                      names(unlist(listdf)))])])), col = k, xaxt = "n")
    points(x=xs, y=listdf[,k+1], cex=0.45)
    legend("topleft", c(names(listdf)[k+1]), 
           col=k, lty=1, cex=0.65)
  }
  #plotting residuals
  resids <- lapply(par$fits, resid)
  plot(y = resids[[k]][1:length(listdf$Cycle)], 
       x = par$fits[[k]]$DATA$Cycles[1:length(listdf$Cycle)], 
       ylim = range(unlist(resids)), 
       xlab = "Cycle", ylab="Fluorescence Residual")
  abline(h=0)
  #finding the RSS
  rss <- sum(resids[[k]]^2)
  #finding the DW-stat
  reg.amp <- dynlm(listdf[,k+1] ~ listdf$Cycle)
  reg.res <- dynlm(resids[[k]] ~ listdf$Cycle)
  dw.amp <- durbinWatsonTest(reg.amp) ; dw.res <- durbinWatsonTest(reg.res)
  #finding the CT value
  ml1 <- modlist(listdf, model = l4)
  res1 <- getPar(ml1, type = "curve", cp = "cpD2", eff = "sliwin")
  #combining all into df
  values <- data.frame(apply(par$params[k,], c(1,2), as.numeric))
  values[,(length(values)+1):(length(values)+9)] <- c(dw.amp$r, dw.amp$dw, dw.amp$p,
                                                      dw.res$r, dw.res$dw, dw.res$p,
                                                      rss, res1[,k][1], res1[,k][2])
  names(values) <- c(names(par$params), paste0(names(dw.amp)[1:3], "-amp"),
                     paste0(names(dw.res)[1:3], "-res"), c("rss", "dw", "p"))
  return(values)
}

subsubplots_tog(l4, try[[5]], 4, plot=T)

subsubplots_tog(l5, try[[5]], 4, plot=T)


##subplot for all
#not nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst)

#nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams(tst)

source("GAPDH.SO/plot_sig.R")
plot_sig(l4, try, plot=TRUE)

plot_sig(b5, try, plot=TRUE)

source("GAPDH.SO/genparams.R")
result.try <- sub_genparams(b5, try[[3]])
KW1_tstresults <- sub_genparams(b5, try[[4]])

source("GAPDH.SO/plot_sig.R")
plot_sig(b5, try, plot=TRUE)
plot_sig(b5, try)
plot_sig(l5, try, macro=1, z=4, plot=TRUE) #tester for NA point
plot_sig(b5, try, macro=3, z=4, plot=TRUE) #tester for good point

####Verification -- creating two matrices 
library(miRcompData) ; data("miRcompData")
targnames <- unique(c(miRcompData$TargetName)) 
miRcompData2 <- miRcompData
miRcompData2 <- miRcompData2[!(miRcompData2$SampleID == "E1_" | miRcompData2$SampleID == "E2_"),]
source("GAPDH.SO/targlist.R")
targetatt <- singtarget.list(miRcompData2, target = targnames[1])

for(i in 1:10){
  for(j in 1:4){
    names[[i]][[j]] <- tst[[i]][[j]]$SampleID
  }
}

strsplit(tst[[1]][[1]]$SampleID[[1]], "_")[[1]][2]

#new attempt
files <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets", 
             pattern =  "^targ_")
for(k in 1:length(files)){
  load(file = files[[k]])
}
targnames <- unique(c(miRcompData2$TargetName)) 
sampnames <- sort(unique(c(miRcompData2$SampleID)))

tmp <- rep(sampnames, 758)
res <- data.frame(
  TargetName = rep(targnames, each = 40), 
  SampleID = rep(sampnames, 758), 
  Group = gsub("_." , "", tmp), 
  FeatureSet = rep(unique(tst[[1]][[1]]$FeatureSet), each = 40), 
  b = rep(NA, 758 * 40), 
  c = rep(NA, 758 * 40), 
  d = rep(NA, 758 * 40), 
  e = rep(NA, 758 * 40), 
  f = rep(NA, 758 * 40), 
  r-amp = rep(NA, 758 * 40), 
  dw-amp = rep(NA, 758 * 40), 
  p-amp = rep(NA, 758 * 40), 
  r-res = rep(NA, 758 * 40), 
  dw-res = rep(NA, 758 * 40), 
  p-res = rep(NA, 758 * 40), 
  rss = rep(NA, 758 * 40), 
  ct = rep(NA, 758 * 40), 
  eff = rep(NA, 758 * 40)
)











