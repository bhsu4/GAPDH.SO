#### Figure 1: Phases of qPCR

library(qpcR) #load in qpcR package for sigmoidal models
source("genparams.R") #generating parameters of sigmoidals

#the plot will use the following target 
load("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targets/targ_hsa-miR-520e_001119.Rda")
miR520e <- unlist.genparams(tst)
nice_b5 <- sub_genparams(est=b5, listdf=miR520e$E) #miR520e set E sigmoidal models, we choose E2

#standardize fluorescence
miR520e.E2 <- miR520e[[5]][,2] 
miR520e.stdE2 <- (miR520e.E2 - min(miR520e.E2))/(max(miR520e.E2) - min(miR520e.E2))

#plot Fluorescence v. Cycle
plot(miR520e.stdE2, xlab = "", ylab='Fluorescence')
mtext("Cycle", side = 1, line = 2)
abline(h=0.045, lty= "dashed")
points(21.5, 0.045, pch=15, cex = 2.0)
segments(21.5, 0.045, 21.5, min(par("usr")))
#arrows(21,0.095,21.35,0.055, length=0.12, angle= 50, lty=3) 
text(23.25, 0.021, "CT Value", col = "black", cex=1)
text(45,0.066, "Threshold", cex=1.05)
mtext('Phase 1: Baseline', side=3, line=-7.5, at=9, font=2)
mtext('Phase 2: Exponential', side=3, line=-7.5, at=24.5, font=2)
mtext('Phase 3: Plateau', side=3, line=-7.5, at=38, font=2)
library(pBrackets)
brackets(26.25, 0.30, 25.25, .21, lwd=1, curvatur=1, type=2)
text(29.75, 0.215, expression(frac('F'[n+1],'F'[n]) == 'PCR Efficiency'), cex=1)
polygon(c(min(par("usr")), 20, 20, min(par("usr")), min(par("usr"))), 
        c(min(par("usr")), min(par("usr")), max(par("usr")), 
          max(par("usr")), min(par("usr"))),
        col= rgb(0,0,0,alpha=0.15), border = NA) 
polygon(c(20,29,29,20,20), 
        c(min(par("usr")), min(par("usr")), max(par("usr")), 
          max(par("usr")), min(par("usr"))),
        col= rgb(0,0,0,alpha=0.12), border=NA) 
polygon(c(29,max(par("usr")),max(par("usr")),29,29 ), 
        c(min(par("usr")), min(par("usr")), max(par("usr")), 
          max(par("usr")), min(par("usr"))),
        col= rgb(0,0,0,alpha=0.09), border= NA)  

text(3.5,0.6125,"Baseline w/ slow upward trend. \n\nObserved fluorescence is usually random\nnoise. The signal is dominated by\nbackground fluorescence.", pos=4)
text(20.5,0.6,"Upward swing into an\nexponential model.\n\nSignal dominated by\nactual amount of\ntarget molecule.", pos=4)
text(32,0.63,"Enters plateau as the signal tapers off.\n\nReaction limited by amount of available\nnucleotides, slowing down of amplification.", pos=4)



#### Figure 2: Residuals for representative curve vs. replicates within subset

setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/GAPDH.SO")
library(qpcR)
source("genparams.R")
source("resid_subplot.R")

#use GAPDH.SO subset and baseline-subtracted subset
setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/")
df <- read.csv(file="GAPDH.SO.csv", header = TRUE, sep = ",")
Cycle = c(1:40)
library(data.table)
setnames(df, "F6.1", "F7") #changing one column name
subsets <- lapply(LETTERS[1:8], function(k) cbind(Cycle,df[,colnames(df) %like% k]))
#colnames(subsets$F) <- c("Cycle", paste0("F", 1:12))
names(subsets) <- LETTERS[1:8]
names(subsets$F)[[8]] <- 'F7'
#list2env(subsets, envir=.GlobalEnv)
subsets$C <- subsets$C[,-c(1)]

#GAPDH.SO data run sigmoidal models
df_b5 <- genparams(est=b5, listdf=subsets)
df_l5 <- genparams(est=l5, listdf=subsets)
df_b4 <- genparams(est=b4, listdf=subsets)
df_l4 <- genparams(est=l4, listdf=subsets)

#baseline subtracted gapdhso
subsets_base <- vector('list', 8) ; subsets2 <- list()
#baseline subtract GAPDH
for(i in 1:8){
  for(j in 2:13){
    subsets_base[[i]][[1]] <- subsets[[i]][,1]
    subsets_base[[i]][[j]] <- subsets[[i]][,j] - min(subsets[[i]][,j])
  }
  subsets2[[i]] <- data.frame(matrix(unlist(subsets_base[[i]]), ncol=13))
}
names(subsets2) <- LETTERS[1:8]
my_names <- lapply(LETTERS[1:8], function(x) paste0(x, 1:12))
for(i in 1:8){
  names(subsets2[[i]]) <- c("Cycle", my_names[[i]])
}

#GAPDH.SO data run sigmoidal models baseline subtracted
df_b5_base <- genparams(est=b5, listdf=subsets2)
df_l5_base <- genparams(est=l5, listdf=subsets2)
df_b4_base <- genparams(est=b4, listdf=subsets2)
df_l4_base <- genparams(est=l4, listdf=subsets2)

dev.off()
m <- matrix(c(1,2,3,4,5,6,7,8),nrow = 2, ncol = 4, byrow = TRUE)
layout(mat = m, widths=c(100, 100, 100, 100))
par(oma=c(4,4,0.5,4),mar=c(0.25,0.25,0,0))
#top 4: gapdh.so
resid_rep_8top(df_l5, mod=l5, 8)
resid_rep_8top(df_b5, mod=b5, 8)
resid_rep_8top(df_l4, mod=l4, 8)
resid_rep_8top(df_b4, mod=b4, 8)
#bottom 4: gapdh.so w/ baseline subtraction
resid_rep_8bot(df_l5_base, l5, 8)
resid_rep_8bot(df_b5_base, b5, 8)
resid_rep_8bot(df_l4_base, l4, 8)
resid_rep_8bot(df_b4_base, b4, 8)
mtext(text="Residual", side=2, line=2, outer=T)
mtext(text="Cycle", side=1, line=2, outer=T)



#### Figure 3: Residuals for each replicate curve vs. actual values

#using same subsets data set (GAPDH.SO) in Figure 2. And calling on same source functions.

#pick one replicate G of GAPDH.SO -- find residuals with each replicate as curve of its own
b5_resultsG <- sub_genparams(b5, subsets$G)
b4_resultsG <- sub_genparams(b4, subsets$G)
l5_resultsG <- sub_genparams(l5, subsets$G)
l4_resultsG <- sub_genparams(l4, subsets$G)

#Run STAR model on GAPDH data set w/ lag = 1, and lag = 2
lstar_mod1 <- list()
for (j in 1:8) lstar_mod1[[j]] <- lapply(subsets[[j]][,2:13], function(f) try(tsDyn::lstar(f, m=1, d=1)))
#lag2
lstar_mod2 <- list()
for (j in 1:8) lstar_mod2[[j]] <- lapply(subsets[[j]][,2:13], function(f) try(tsDyn::lstar(f, m=2, d=1)))

#format of plot
dev.off()
m <- matrix(c(1,2,5,3,4,6),nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, widths=c(100, 100, 100))
par(oma=c(4,4,0.5,4),mar=c(0.25,0.25,0,0))

#call on subplots
subplot_resid_5xy(l5_resultsG)
subplot_resid_5x(b5_resultsG)
subplot_resid_4xy(l4_resultsG)
subplot_resid_4x(b4_resultsG)

for(i in 1:12){
  if(i==1){
    plot(x=2:40, y=lstar_mod1[[7]][[i]]$residuals, type = 'l', col = 1, 
         ylim = c(-11000,9000), xaxt = "n", yaxt = "n")
    axis(side=4, at=seq(-8000, 8000, by=4000)) }
  else{ lines(x=2:40, y=lstar_mod1[[7]][[i]]$residuals, col=i) }
}

for(i in 1:12){
  if(i==1){
    plot(x=3:40, y=lstar_mod2[[7]][[i]]$residuals, type = 'l', 
         col = 1, ylim = c(-11000, 9000), xaxt = "n", yaxt = "n")
    axis(side=4, at=seq(-8000, 8000, by=4000))
    axis(side=1, at=seq(0, 40, by=10))
  }
  else{ lines(x=3:40, y=lstar_mod2[[7]][[i]]$residuals, col=i) }
}
mtext(text="Cycle", side=1, line=2, outer=TRUE)
mtext(text= "Residual", side=2, line=2, outer=TRUE)



#### Figure 4: DW Test Statistic for Residuals

library(ggplot2)
library(miRcompData)
data("miRcompData")
miRcompData2 <- miRcompData #we will use the following large data set
miRcompData2 <- miRcompData2[!(miRcompData2$SampleID == "E1_" | miRcompData2$SampleID == "E2_"),] #delete E's

#loading in matrices of parameter est
setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targetsmcont")
getmat <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targetsmcont", 
              pattern = "_dat")

#TO OBTAIN: the sigmoidal parameters and statistics -- use sig_est
source("GAPDH.SO/plot_sig.R")

#larger data set: with miRcomp
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont", 
                pattern =  "^targ_")

#running to find matrix combined
l5dat2 <- sig_est(l5, miRcompData2, getfiles=getfiles) #this will take a substantial amount of time
save(l5dat2, file = paste0("l5_" , "dat2", ".Rda")) #save file

#load our file in
load(file = getmat[[8]])

#heat map in contour
par(oma=c(2,2,0.5,0.5),mar=c(1.25,2,0.75,0.75),mfrow=c(2,2),pch=16)
commonTheme = list(labs(color="Density",fill="Density",
                        x="Durbin-Watson Test Statistic",
                        y="Cq Value (SDM)"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))
#dw.comp
l5dat2.na <- l5dat2[which(  (!is.na(l5dat2$dw.comp)) & (!is.na(l5dat2$ct))),]

ggplot(data=l5dat2.na,aes(dw.comp, ct)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  # geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  xlim(c(0,2.5))+ #ylim(c(5,40)) +
  guides(alpha="none") +
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme



#### Figure 5: Classification of Different Signals

setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/GAPDH.SO")
library(dynlm) ; library(car) ; library(qpcR) ; library(Hmisc)
source("genparams.R")
source("plot_sig.R")
source("resid_subplot.R")
source("l5_model.R") ; source("b5_model.R") ; source("l4_model.R") ; source("b4_model.R")

##good, bad, non-signals, unidentified signals (failed exp)
#nice dataset
setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targets")
load("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams(tst) #nice clean data set
load("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst) #not nice data set with more noise

#run this to see all curves: for specified set of data -- we will apply this same function for our figure
plot_sig(l5, try, plot=TRUE) #plots shown
plot_sig(l5, try.good, plot=TRUE) #plots shown
plot_sig(l4, try.good, plot=TRUE) #plots shown

#how we generated our figure: plot_sig is general application, sigmacro is an adhoc variation of plot_sig for the figure
par(oma=c(2,5,0,3), mar=c(0,0,1,2), mfrow=c(5,2),pch=16)
plot_sigmacro(l5, try.good, macro=4, z=8, plot=TRUE) #goodfit nicedata l5
plot_sigmacro(l4, try.good, macro=4, z=9, plot=TRUE) #poorfit nicedata l4
plot_sigmacro(l5, try.good, macro=4, z=9, plot=TRUE) #poorfit nicedata l5
plot_sigmacro(l5, try, macro=4, z=3, plot=TRUE) #nonsignal notnicedata l5
plot_sigmacro(l5, try, macro=2, z=9, plot=TRUE) #failed notnicedata l5
mtext(text="Cycle", side=1, line=1, outer=T)
mtext(text="Fluorescence", side=2, line=3, outer=T)
mtext(text="Residual", side=4, line=1, outer=T)



#### Table 1: Residual Sum of Squares and Local RSS for L5 and LSTAR

#larger data set: with miRcomp
#loading in matrices of parameter est
setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targetsmcont")
getmat <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targetsmcont", 
              pattern = "_dat")
load(file = getmat[[9]]) #l5
getmat <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targetsmcont", 
              pattern = "LSTAR")
load(file = getmat[[3]]) #lstar lag2

#our df of results are l5dat.2, and tstlstarmat3 -- rss (Residual Sum of Squares), rssred (Local RSS)
summary(l5dat.2[,c("rss", "rssred")])
summary(tstlstarmat3[,c("rss", "rssred")])



#### Supplementary Figure 1: AR(1), AR(2) Models

load("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targets/targ_hsa-miR-520e_001119.Rda")
miR520e <- unlist.genparams(tst)
miR520e.E2 <- miR520e[[5]][,2] #use E2

par(mfrow = c(1,2))
par(oma=c(4,4,0.5,0.5),mar=c(0.25,0.25,0,0))

#ar(1) model 
artest <- arima(x = ts(miR520e.E2), order = c(1,0,0), method = 'ML')
plot(miR520e.E2, xlab = "Cycle", ylab="Fluorescence")
lines(2:46, (miR520e.E2*artest$coef[1])[1:45], col=1)
#ar(2)
artest2 <- arima(x = ts(miR520e.E2), order = c(2,0,0), method = 'ML')
plot(miR520e.E2, xlab = "Cycle", ylab="", yaxt = "n")
lines(3:46, (miR520e.E2*artest2$coef[1])[1:44]+(miR520e.E2*artest2$coef[2])[1:44], col=1)
mtext("Cycle", line = 2, side = 1, outer = TRUE)
mtext("Fluorescence", line = 2.5, side = 2, outer = TRUE)



#### Supplementary Figure 2: AIC to choose Lag 2

library(dynlm) ; library(Hmisc)
setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/GAPDH.SO")
source("lag_gen.R")

#time-series conversion
for(i in 1:8){
  for(j in 1:13){
    ts_subsets[[i]][,j] <- ts(subsets[[i]][,j])
  }
}
fftest <- lapply(LETTERS[1:8], paste0, 1:12) #list of subset names
fftry <- lapply(fftest, function(x) lapply(x, get_lag_formulae, n=15)) #dynlm formula in list
ff_fitstest <- list()
for(j in 1:8) ff_fitstest[[j]] <- lapply(fftry[[j]], function(x) lapply(x, 
                                                     function(f) dynlm(formula = as.formula(f), 
                                                                          data=ts_subsets[[j]]))) #output results of dynlm
ff_fitstestaic <- list()
for(j in 1:8) ff_fitstestaic[[j]] <- sapply(ff_fitstest[[j]], function(x) sapply(x,AIC)) #output results of dynlm

#plot our AICs
dev.off()
for(j in 1:8){
  for(i in 1:12){
    if(i==1 & j==1){
      plot(x=1:15, y=ff_fitstestaic[[j]][,i], type = "p", ylim = c(range(ff_fitstestaic)), 
           ylab = "AIC", xlab = "Lag Term")
    }
    points(x=1:15, y=ff_fitstestaic[[j]][,i], col = j)
  } #plotting all replication sets' AIC
}
for(j in 1:8){
  lines(spline(x = 1:15, y=rowMeans(ff_fitstestaic[[j]])), pch = 19, cex = 3, col = j)
} #smoothing line



#### Supplementary 3: miRcomp residuals between replicate residuals

setwd("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/GAPDH.SO")
source("resid_subplot.R")
#load in target
load("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")

#Raw non baseline-subtracted
miR23a.Rn <- unlist.genparams.Rn(tst)
miR23a.Rn_b5 <- genparams(est=b5, listdf=miR23a.Rn)
miR23a.Rn_l5 <- genparams(est=l5, listdf=miR23a.Rn)
miR23a.Rn_b4 <- genparams(est=b4, listdf=miR23a.Rn)
miR23a.Rn_l4 <- genparams(est=l4, listdf=miR23a.Rn)

#baseline-subtracted
miR23a.dRn <- unlist.genparams(tst)
miR23a.dRn_b5 <- genparams(est=b5, listdf=miR23a.dRn)
miR23a.dRn_l5 <- genparams(est=l5, listdf=miR23a.dRn)
miR23a.dRn_b4 <- genparams(est=b4, listdf=miR23a.dRn)
miR23a.dRn_l4 <- genparams(est=l4, listdf=miR23a.dRn)

#format plot
dev.off()
m <- matrix(c(1,2,3,4,5,6,7,8),nrow = 2, ncol = 4, byrow = TRUE)
layout(mat = m)
par(oma=c(4,4,4,4),mar=c(0,0.5,0.75,0))

#residuals of raw fluoresence
resid_rep_8top_miR(miR23a.Rn_l5, mod=l5, 6)
resid_rep_8top_miR(miR23a.Rn_b5, mod=b5, 6)
resid_rep_8top_miR(miR23a.Rn_l4, mod=l4, 6)
resid_rep_8top_miR(miR23a.Rn_b4, mod=b4, 6)
#residuals of baseline-subtr fluoresence
resid_rep_8bot_miR(miR23a.dRn_l5, l5, 6)
resid_rep_8bot_miR(miR23a.dRn_b5, b5, 6)
resid_rep_8bot_miR(miR23a.dRn_l4, l4, 6)
resid_rep_8bot_miR(miR23a.dRn_b4, b4, 6)
mtext(text="Residual", side=2, line=2, outer=T)
mtext(text="Cycle", side=1, line=2, outer=T)



#### Supplementary Figure 4: miRcomp residuals for each replicate vs. actual 

source("resid_subplot.R")
#load target
load("C:/Users/Benjamin Hsu/Desktop/Other/Fall 2018/Independent Study/GAPDH.SO/targets/targ_hsa-miR-500_002428.Rda")
try <- unlist.genparams(tst)

#lag1
lstar_miR1 <- list()
for (j in 1:10) lstar_miR1[[j]] <- lapply(try[[j]][,2:5], function(f) try(tsDyn::lstar(f, m=1, d=1)))
#lag2
lstar_miR2 <- list()
for (j in 1:10) lstar_miR2[[j]] <- lapply(try[[j]][,2:5], function(f) try(tsDyn::lstar(f, m=2, d=1)))

#format plots
dev.off()
m <- matrix(c(1,2,5,3,4,6),nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, widths=c(100, 100, 100))
par(oma=c(4,4,0.5,4),mar=c(0.25,0.25,0,0))

#plot
subplot_resid_miR(try$H) #sigmoidal plots 
#lag1
for(i in 1:4){
  if(i==1){
    plot(x=2:46, y=lstar_miR1[[8]][[i]]$residuals, type = 'l', 
         col = 1, ylim = c(-150,150), #c(range(lapply(lstarmodnolog1[[7]], function(x) unlist(x$residuals)))),
         xaxt = "n", yaxt = "n")
    axis(side=4, at=seq(-100, 100, by=50)) }
  else{ lines(x=2:46, y=lstar_miR1[[8]][[i]]$residuals, col=i) }
}
#lag2
for(i in 1:4){
  if(i==1){
    plot(x=3:46, y=lstar_miR2[[8]][[i]]$residuals, type = 'l', 
         col = 1, ylim = c(-150, 150),   #range(lstarmodnolog2[[7]], na.rm=TRUE)), 
         xaxt = "n", yaxt = "n")
    axis(side=4, at=seq(-100, 100, by=50))
    axis(side=1, at=seq(0, 40, by=10))  }
  else{ lines(x=3:46, y=lstar_miR2[[8]][[i]]$residuals, col=i) }
}



#### Supplementary Figure 5: Durbin-Watsson, Ljung-Box and Pearson Correlation

#load our file in: Same file as in Figure4
load(file = getmat[[8]])
l5dat2.na <- l5dat2[which(  (!is.na(l5dat2$dw.comp)) & (!is.na(l5dat2$ct))),]

#ljung-box p-value vs. durbin-watson
commonTheme = list(labs(color="Density",fill="Density",
                        x="Durbin-Watson Test Statistic",
                        y="Ljung-Box P-value"),
                   theme_bw(), theme(legend.position = "none"))
plot3 <- ggplot(data=l5dat2.na,aes(dw.comp,boxlj.p)) + 
  ylim(0,1) + xlim(c(0,4)) + guides(alpha="none") + 
  geom_hline(yintercept=0.05, lty='dashed') +
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme

#pearson corr vs. ljung-box p-value
commonTheme = list(labs(color="Density",fill="Density",
                        x="Ljung-Box P-value",
                        y="Pearson Correlation"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))
plot4 <- ggplot(data=l5dat2.na,aes(boxlj.p,pcor)) + 
  #stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  #scale_fill_continuous(low="green",high="red") +
  #geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  ylim(-1,1)+
  xlim(c(0,1))+
  guides(alpha="none") + geom_vline(xintercept=0.05, lty='dashed')+
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme

grid.arrange(plot3, plot4, ncol=2)

