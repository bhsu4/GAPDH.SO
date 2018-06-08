plot_resid <- function(listdf, params) {
  
  for (i in 1:length(params$fits)){
    for(k in 1:(length(listdf[[i]])-1)){
      resids <- lapply(params$fits, resid)
      if(k == 1) plot(y=resids[[i]][1:length(listdf[[i]]$Cycle)], 
                      x=params$fits[[i]]$DATA$Cycles[1:length(listdf[[i]]$Cycle)], 
                      ylim=range(resids[[i]]), type="l", 
                      xlab="Cycle", ylab="Fluoresence Residual")
      if(k > 1){
        ind2 <- length(listdf[[i]]$Cycle)*k
        ind1 <- ind2-(length(listdf[[i]]$Cycle)-1)
        lines(y=resids[[i]][ind1:ind2], x=params$fits[[i]]$DATA$Cycles[ind1:ind2], col=k)
      }
    }
    title(main= paste(names(params$fits[i]), sub=params$fits$A$MODEL$name, sep = ", "))
  }
}

subplot_resid <- function(params){
  
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1) plot(y=resids[[k]][1:40], 
                    x=params$fits[[1]]$DATA$Cycles[1:40], 
                    ylim=range(resids[[k]]), type="l", 
                    xlab="Cycle", ylab="Fluorescence")
    if(k > 1){
      ind2 <- 40*k
      ind1 <- ind2-39*(k)-(k-1)
      lines(y=resids[[k]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
}
