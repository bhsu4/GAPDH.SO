setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/")
df <- read.csv(file="GAPDH.SO.csv", header = TRUE, sep = ",")
Cycle = c(1:40)
library(data.table)
subsets <- lapply(LETTERS[1:8], function(k) cbind(Cycle,df[,colnames(df) %like% k]))
names(subsets) <- LETTERS[1:8]
#list2env(subsets, envir=.GlobalEnv)
subsets$C <- subsets$C[,-c(1)]

library(qpcR)
source("GAPDH.SO/genparams.R")
df_b5 <- genparams(est=b5, listdf=subsets)
df_l5 <- genparams(est=l5, listdf=subsets)
df_b4 <- genparams(est=b4, listdf=subsets)
df_l4 <- genparams(est=l4, listdf=subsets)

## to get residuals
source("GAPDH.SO/plot_resid.R")
b5resids <- lapply(df_b5$fits, resid)
plot_resid(df_b5, b5resids)

source("GAPDH.SO/b5_model.R")
plot_b5(Cycle, subsets, df_b5)
source("GAPDH.SO/b4_model.R")
plot_b4(Cycle, subsets, df_b4)
source("GAPDH.SO/l4_model.R")
plot_l4(Cycle, subsets, df_l4)
source("GAPDH.SO/l5_model.R")
plot_l5(Cycle, subsets, df_l5)

source("GAPDH.SO/plot_subset.R") #by model
plot_subset(Cycle, df_l4, df_l5, df_b4, df_b5, subsets)
source("GAPDH.SO/plot_together.R") 
plot_together(Cycle, df_l4, df_l5, df_b4, df_b5, subsets)

### subsetting with efficiencies / run fluo eff plots

source("GAPDH.SO/eff_phase.R")
Cycles = 1:39
subsetsb_b5 <- genparamsbase(subsets, df_b5, subsets)

#creating ground phase slanted baseline model
plot_eff(Cycles, subsetsb_b5, subsets, subsets, "slant", 11)
plot_eff(Cycles, subsetsb_b5, subsets, subsets, "flat", 11)
plot_eff(Cycles, subsetsb_b5, subsets, curve_b5, "slant", 11)

#### confirmation on curve values fitted

testconf <- function(x, b, e, f){
  ((1+exp(b*((x-1)-e)))/(1+exp(b*(x-e))))^f
}

plot_testconf <- function(xs, listdf, par){
  plot(x=xs, y=testconf(xs, b=par$params$b[1], e=par$params$e[1], f=par$params$f[1]), type="l",  
       xlab="Cycle", ylab="Fluorescence")
  for(k in 2:length(subsets)){
    lines(x=xs, y=testconf(xs, b=par$params$b[k], e=par$params$e[k], f=par$params$f[k]), col=k)
  }
  legend("topright", c(LETTERS[1:length(subsets)]), 
         col=1:length(subsets), ncol=2, lty=1, cex=0.5)
}

plot_testconf(Cycles, subsets, df_b5) #looks like eff. from curve values

##Applying STAR Model

#Chow test single break confirm
source("GAPDH.SO/strucbreak.R")
breakA <- strbreak_chow(Cycle, subsets)

plot(x=5:35, y = breakA[[1]][[1]][,1], ylab = "F-statistic", xlab = "Cycle")
for(k in 1:8){
  for(i in 1:12){
    lines(x=5:35, y=breakA[[1]][[k]][,i], col=k)
  } #breakpoints for replicates
}
for(j in 2:ncol(subsets$A)){
  for(h in 1:length(subsets)){
    points(x=Cycle, y=(50*subsets[[h]][,j]/max(subsets[[h]][,j])), cex=0.45)
  }
} #scaled overlay of points

#creating list for values to run through nrep
names(ts_subsets$F)[8] <- "F7" #replace F.6.1

fftest <- lapply(LETTERS[1:8], paste0, 1:12) #list of subset names
fftry <- lapply(fftest, function(x) lapply(x, get_lag_formulae, n=15)) #dynlm formula in list
ff_fitstest <- list()
for(j in 1:8) ff_fitstest[[j]] <- lapply(fftry[[j]], function(x) lapply(x, 
                                                     function(f) dynlm(formula = as.formula(f), 
                                                     data=ts_subsets[[j]]))) #output results of dynlm
for(j in 1:8) ff_fitstestaic[[j]] <- sapply(ff_fitstest[[j]], 
                                            function(x) sapply(x,AIC)) #output results of dynlm

#n-lag AIC values
for(j in 1:8){
  for(i in 2:16){
    if(j==1){
      plot(x=1:12, y=ff_fitstestaic[[j]][1,], type = "p", 
           ylim=range(ff_fitstestaic), col = i, xlab = "Replication Set", ylab = "AIC")
    }
    points(x=1:12, y=ff_fitstestaic[[j]][i-1,], col = i-1)
  }
}

#finding 5% and 95% plot/ deciding on lag term

for(j in 1:8){
  for(i in 1:12){
    if(i==1 & j==1){
      plot(x=1:15, y=ff_fitstestaic[[j]][,i], type = "p", ylim = c(range(ff_fitstestaic)), 
           ylab = "AIC", xlab = "Replication Set")
    }
    points(x=1:15, y=ff_fitstestaic[[j]][,i], col = j)
  } #plotting all replication sets' AIC
}
for(j in 1:8){
  lines(spline(x = 1:15, y=rowMeans(ff_fitstestaic[[j]])), pch = 19, cex = 3, col = j)
} #smoothing line
qts <- quantile(range(ff_fitstestaic), probs=c(0.05, 0.95, 0.15, 0.85)) #percentage lines of range
abline(h=qts[1], col="red") 
abline(h=qts[2], col="red")
abline(h=qts[3], col="green")
abline(h=qts[4], col="green")
#use lag(2) since dramatically better, afterwards plateaus

##Bai and Perron's Test
A1test <- log10(subsets$A$A1)
A1test <- cbind(A1test, lag(subsets$A$A1, k=-2))
colnames(A1test) <- c("y", "ylag2")
A1test <- ts(A1test)
re.a1 <- efp(y ~ ylag1, data = A1test, type = "RE")
plot(re.a1)
bp.a1 <- breakpoints(y ~ ylag2, data = A1test)
summary(bp.a1)
plot(A1test[,"y"], ylab = expression(log[10](Fluorescence))) 
lines(bp.a1, breaks = 2)

#breakpoints at each replication
subsets_log <- subsets
bp.subs <- vector("list", 12)
for(i in 1:8){
  for(j in 2:13){
    subsets_log[[i]][[j-1]] <- log10(subsets[[i]][[j]]) #converting to log10
    subsets_log[[i]][[j-1]] <- cbind(subsets_log[[i]][[j-1]], lag(subsets[[i]][[j]], k=-2)) #include lag2
    colnames(subsets_log[[i]][[j-1]]) <- c("y", "ylag2")
    subsets_log[[i]][[j-1]] <- ts(subsets_log[[i]][[j-1]]) #change to time series to apply brkpt
    bp.subs[[i]][[j-1]] <- breakpoints(y ~ ylag2, data=subsets_log[[i]][[j-1]]) #finding brkpt for each
    if(j==2){
      plot(subsets_log[[i]][[j-1]][,"y"], ylab = expression(log[10](Fluorescence)), 
           ylim = c(range(subsets_log[[i]][[j-1]][,"y"])), col = i) 
    } #plotting each subset A-H
    lines(subsets_log[[i]][[j-1]][,"y"], ylab = expression(log[10](Fluorescence)),
          ylim = c(range(subsets_log[[i]][[j-1]][,"y"])), col = i)
    lines(bp.subs[[i]][[j-1]], breaks = 2, col = j-1) #adding brkspt cycle times
  }
}

##plotting LSTAR model to see fit for replication
df <- data.frame(1:40)
ff <- rowMeans(lsubsets$A[,2:13])
ff2 <- apply(lsubsets$A[,2:13], 1, median)

try.lstar <- tsDyn::lstar(ff, m=2, d=1) #embedding dimension=2, delay = 1 #mean
plot(x=1:38, try.lstar$fitted.values, type = "l", ylim = c(5.1, 5.75))
for(j in 2:13){
  for(h in 1:8){
    points(x=1:38, y=subsets_log[[h]][[j]][1:38,1], cex=0.45)
  }
} #plot points for 38 cycles given lag 2
try.lstar2 <- tsDyn::lstar(ff2, m=2, d=1) #median replication set
lines(x=1:38, try.lstar2$fitted.values, type = "l", ylim = c(5.1, 5.75), col = 2)

df_b5_log <- genparams(est=b5, listdf=lsubsets)
lines(x=1:38, y=b5_model(1:38, b=df_b5_log$params$b[1], c=df_b5_log$params$c[1],
                         d=df_b5_log$params$d[1], e=df_b5_log$params$e[1], 
                         f=df_b5_log$params$f[1]), col=2) #lines for b_5 model

