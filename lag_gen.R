power <- function(n){
  function(x) x^n
}



get_lag_formulae <- function(n, nrep){
  formulae <- list()
  for(k in 1:n){
    #A1Test <- ts(subsets$A$A1) run list equivalence, then run subsets thru it creating new
    formulae[[k]] <- paste(nrep, "~", paste(paste0("L(", nrep, ",", 1:k,")"), collapse=" + "), collapse=" ")
  }
  return(formulae)
}

ff <- get_lag_formulae(6, nrep="A5")

ff_fits <- lapply(ff, function(f) dynlm(formula = as.formula(f), data=subsets[[1]]))
sapply(ff_fits, AIC)


anova.dyn(dynlm(subsets$A$A1 ~ L(subsets$A$A1)))

data("Nile")
dynlm(Nile ~ L(Nile, 1))
dynlm(Nile ~ L(Nile, 1:3))
dynlm(A1Test ~ L(A1Test, 1))
