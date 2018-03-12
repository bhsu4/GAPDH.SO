plot_resid <- function(params,resids){

for (i in 1:length(params$fits)){
  for(k in 1:12){
    resids <- lapply(params$fits, resid)

    if(k == 1) plot(y=resids[[i]][1:40], 
                    x=params$fits[[1]]$DATA$Cycles[1:40], 
                    ylim=range(resids[[i]]), type="l", 
                    xlab="Cycle", ylab="Fluoresence Residual")
    if(k > 1){
      ind2 <- 40*k
      ind1 <- ind2-39
      lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
title(main= names(params$fits[i]))
  }
}