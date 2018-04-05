power <- function(n){
  function(x) x^n
}

#A1Test <- ts(subsets$A$A1) run list equivalence, then run subsets thru it creating new

get_lag_formulae <- function(n, nrep){
    ts_subsets <- subsets
  for(i in 1:8){
  for(k in 2:13){
    ts_subsets[[i]][[k]] <- ts(subsets[[i]][[k]])
    }
  }
    formulae <- list()
  for(k in 1:n){
    formulae[[k]] <- paste(nrep, "~", paste(paste0("L(", nrep, ",", 1:k,")"), collapse=" + "), collapse=" ")
  }
  return(formulae)
}

#creating list for values to run through nrep
names(ts_subsets$F)[8] <- "F7"

fftest <- lapply(LETTERS[1:8], paste0, 1:12) #list of subset names
fftry <- lapply(fftest, function(x) lapply(x, get_lag_formulae, n=15)) #dynlm formula in list
ff_fitstest <- list()
for(j in 1:8) ff_fitstest[[j]] <- lapply(fftry[[j]], function(x) lapply(x, 
                                function(f) dynlm(formula = as.formula(f), 
                                                  data=ts_subsets[[j]]))) #output results of dynlm
for(j in 1:8) ff_fitstestaic[[j]] <- sapply(ff_fitstest[[j]], function(x) sapply(x,AIC)) #output results of dynlm

#looking at 4 lagged AIC values
for(j in 1:8){
  for(i in 2:5){
if(i==1){
plot(x=1:12, y=ff_fitstestaic[[j]][1,], type = "l", 
     ylim=range(ff_fitstestaic), col = i, xlab = "Lagged Parameter", ylab = "AIC")
}
    lines(x=1:12, y=ff_fitstestaic[[j]][i-1,], col = i-1)
  }
}

#finding 5% and 95% plot

for(j in 1:8){
for(i in 1:12){
  if(i==1 & j==1){
plot(x=1:15, y=ff_fitstestaic[[j]][,i], type = "p", ylim = c(range(ff_fitstestaic)))
  }
points(x=1:15, y=ff_fitstestaic[[j]][,i], col = j)
}
}
for(j in 1:8){
    lines(x = 1:15, y=rowMeans(ff_fitstestaic[[j]]), pch = 19, cex = 3, col = j)
}
qts <- quantile(range(ff_fitstestaic), probs=c(0.05, 0.95, 0.15, 0.85))
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
      }
    lines(subsets_log[[i]][[j-1]][,"y"], ylab = expression(log[10](Fluorescence)),
                                         ylim = c(range(subsets_log[[i]][[j-1]][,"y"])), col = i)
    lines(bp.subs[[i]][[j-1]], breaks = 2, col = j-1)
  }
}



##########finding delay parameter w/ LSTAR

library(tsDyn)

#fit a LSTAR model. Note 'maxit': slow convergence
mod.lstar <- tsDyn::lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
mod.lstar

#fit a LSTAR model without a constant in both regimes. 
mod.lstar2 <- lstar(log10(lynx), m=1,  include="none")
mod.lstar2

#Note in example below that the initial grid search seems to be to narrow. 
# Extend it, and evaluate more values (slow!):
controls <- list(gammaInt=c(1,2000), nGamma=50)
mod.lstar3 <- lstar(log10(lynx), m=1,  include="none", starting.control=controls)
mod.lstar3

# a few methods for lstar:
summary(mod.lstar)
residuals(mod.lstar)
AIC(mod.lstar)
BIC(mod.lstar)
plot(mod.lstar)
predict(mod.lstar, n.ahead=5)
