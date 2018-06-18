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