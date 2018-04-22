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
lsubsets <- subsets #lsubsets different from subsets_log b/c subsets_log with ts 
for(i in 1:8){
  for(j in 2:13){
    lsubsets[[i]][[j-1]] <- log10(subsets[[i]][[j]])
  }
}

ff <- rowMeans(lsubsets$A[,2:13])
ff2 <- apply(lsubsets$A[,2:13], 1, median)

try.lstar2 <- tsDyn::lstar(fftest, m=2, d=1)
try.lstar <- tsDyn::lstar(ff, m=2, d=1) #embedding dimension=2, delay = 1 #mean

plot(x=3:40, try.lstar$fitted.values, type = "l", ylim = c(5.1, 5.75), xlim=c(1,40))
for(j in 1:12){
  for(h in 1:8){
    points(x=1:40, y=subsets_log[[h]][[j]][,1], cex=0.45)
  }
} #plot points for 38 cycles given lag 2try.lstar2 <- tsDyn::lstar(ff2, m=2, d=1) #median replication set
#lines(x=1:38, try.lstar2$fitted.values, type = "l", ylim = c(5.1, 5.75), col = 2)

df_b5_log <- genparams(est=b5, listdf=lsubsets)
lines(x=1:40, y=b5_model(1:40, b=df_b5_log$params$b[1], c=df_b5_log$params$c[1],
                         d=df_b5_log$params$d[1], e=df_b5_log$params$e[1], 
                         f=df_b5_log$params$f[1]), col=2) #lines for b_5 model


#plotting residuals: do we need to forecast, if by looking, resid worse for lstar than log models?
plot(x=1:38, try.lstar$residuals, type = "l") #reduces branching effect?? doesnt match up

ff <- vector("list", 8)
for(i in 1:8){
  ff[[i]] <- rowMeans(lsubsets[[i]][, 2:13])
} #list of 8 for row means

try.lstar <- vector("list", 8)
for(i in 1:8){
  try.lstar[[i]] <- tsDyn::lstar(ff[[i]], m=2, d=1) #run lstar model through each rep
}

for(i in 1:8){
  if(i == 1){ #plot residuals/ note: in ln units
    plot(x = 1:38, y = try.lstar[[i]]$residuals, type = "l", col = 1, cex = 1.5, 
         ylab = "Residuals", xlab = "Time Series Cycle (t)")
  }
  lines(x = 1:38, y = try.lstar[[i]]$residuals, col = i, cex = 1.5)
}
legend("bottomright", c(LETTERS[1:8]), col=1:8, ncol = 4, lty = 1)

#plotting ln units for plot_resid comp

library(qpcR)
source("GAPDH.SO/genparams.R")
df_b5_log <- genparams(est=b5, listdf=lsubsets) #log b5 params
df_b4_log <- genparams(est=b4, listdf=lsubsets) #log b4 params
df_l5_log <- genparams(est=l5, listdf=lsubsets) #log b4 params
df_l4_log <- genparams(est=l4, listdf=lsubsets) #log b4 params

## to get residuals
source("GAPDH.SO/plot_resid.R")
b5lresids <- lapply(df_b5_log$fits, resid)
b4lresids <- lapply(df_b4_log$fits, resid)
l5lresids <- lapply(df_l5_log$fits, resid)
l4lresids <- lapply(df_l4_log$fits, resid)

plot_resid(df_b5_log, b5lresids)
plot_resid(df_b4_log, b4resids)
plot_resid(df_l5_log, l5resids)
plot_resid(df_l4_log, l4resids)


##each replication set 
lstarmod <- list()
for(j in 1:8) lstarmod[[j]] <- lapply(lsubsets[[j]][,2:13], function(f) try(tsDyn::lstar(f, m=2, d=1)))
#doesn't work???

lstarmod1 <- lapply(lsubsets[[1]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))
lstarmod2 <- lapply(lsubsets[[2]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))
lstarmod3 <- lapply(lsubsets[[3]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))
lstarmod4 <- lapply(lsubsets[[4]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))
lstarmod5 <- lapply(lsubsets[[5]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))
lstarmod6 <- lapply(lsubsets[[6]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))
lstarmod7 <- lapply(lsubsets[[7]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))
lstarmod8 <- lapply(lsubsets[[8]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))

lstarmod1nolog <- lapply(subsets[[1]][,2:13], function(x) tsDyn::lstar(x, m=2, d=1))

for(h in 1:8){
  for(j in 1:12){
    if(j == 1){
      plot(x=3:40, lstarmod1[[1]]$fitted.values, type = "l", ylim = c(5.1, 5.75), xlim=c(1,40))
    }
    lines(x=3:40, lstarmod1[[j]]$fitted.values, type = "l", col = j)
    
    points(x=1:40, y=subsets_log[[h]][[j]][,1], cex=0.45)
  }
} 

for(i in 1:8){
  if(i == 1){ #plot residuals/ note: in ln units
    plot(x = 1:38, y = lstarmod1[[i]]$residuals, type = "l", col = 1, cex = 1.5, 
         ylab = "Residuals", xlab = "Time Series Cycle (t)")
  }
  lines(x = 1:38, y = lstarmod1[[i]]$residuals, col = i, cex = 1.5)
}
legend("bottomright", c(names(lstarmod1)), col=1:8, ncol = 6, lty = 1, cex=0.75)


#nolog residual plots
for(i in 1:8){
  if(i == 1){ #plot residuals/ note: in ln units
    plot(x = 1:38, y = lstarmod1nolog[[i]]$residuals, type = "l", col = 1, cex = 1.5, 
         ylab = "Residuals", xlab = "Time Series Cycle (t)")
  }
  lines(x = 1:38, y = lstarmod1nolog[[i]]$residuals, col = i, cex = 1.5)
}
legend("bottomright", c(names(lstarmod1nolog)), col=1:8, ncol = 6, lty = 1, cex=0.75)

#original b5
plot_resid(df_b5, b5resids)
