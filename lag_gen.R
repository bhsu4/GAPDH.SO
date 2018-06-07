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
