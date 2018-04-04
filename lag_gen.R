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

fftest <- vector('list', 8)
for(j in 1:8){
for(i in 1:12){
fftest[[j]][[i]] <- paste0(LETTERS[j], i)
  }
}

fftry <- vector('list', 8)
for(j in 1:8){
  for(i in 1:12){
    fftry[[j]][[i]] <- get_lag_formulae(4, fftest[[j]][[i]])
  }
}

ff_fitstest <- vector('list', 8)
for(j in 1:8){
  for(i in 1:12){
ff_fitstest[[j]][[i]] <- lapply(fftry[[j]][[i]], 
                                function(f) dynlm(formula = as.formula(f), data=ts_subsets[[j]]))
  }
}

ff_fitstestaic <- vector('list', 8)
for(j in 1:8){
  for(i in 1:12){
ff_fitstestaic[[j]][[i]] <- sapply(ff_fitstest[[j]][[i]], AIC)
  }
}

plot(x=1, y=ff_fitstestaic[[1]][[1]])
stripchart(ff_fitstestaic, vertical = TRUE)
#plotting x from 12*8 , with y values taking on AIC

#works down here
ff <- get_lag_formulae(4, "A2")
ff_fits <- lapply(ff, function(f) dynlm(formula = as.formula(f), data=ts_subsets[[1]]))
sapply(ff_fits, AIC)
#if only look at 4 lags, then we get use of AR(4)


A1test <- log10(subsets$A$A1)
A1test <- cbind(A1test, lag(subsets$A$A1, k=-4))
colnames(A1test) <- c("y", "ylag4")
A1test <- ts(A1test)
re.a1 <- efp(y ~ ylag4, data = A1test, type = "RE")
plot(re.a1)
bp.a1 <- breakpoints(y ~ ylag4, data = A1test)
summary(bp.a1)
plot(A1test[,"y"], ylab = expression(log[10](Fluorescence))) 
lines(bp.a1, breaks = 2)

mod.lstar <- lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
mod.lstar



##########finding delay parameter w/ LSTAR

library(tsDyn)

#fit a LSTAR model. Note 'maxit': slow convergence
mod.lstar <- lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
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
