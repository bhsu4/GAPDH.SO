subsubplots_tog <- function(est, listdf, k, plot=FALSE){
  #genearting parameters of specific
  par <- try(sub_genparams(est, listdf))
  #two graphs on top each other
  par(mfrow=c(2,1))
  par(oma=c(4,4,4,4),mar=c(0.25,0.25,0,0))
  #finding the residuals
  resids <- lapply(par$fits, resid)
  #finding the RSS
  rss <- sum(resids[[k]]^2)
  #finding the DW-stat
  reg.amp <- dynlm(listdf[,k+1] ~ listdf$Cycle)
  reg.res <- dynlm(resids[[k]] ~ listdf$Cycle)
  dw.amp <- durbinWatsonTest(reg.amp) ; dw.res <- durbinWatsonTest(reg.res)
  #finding the CT value
  ml1 <- modlist(listdf, model = est)
  res1 <- getPar(ml1, type = "curve", cp = "cpD2", eff = "sliwin")
  #combining all into df
  values <- data.frame(apply(par$params[k,], c(1,2), as.numeric))
  values[,(length(values)+1):(length(values)+9)] <- c(dw.amp$r, dw.amp$dw, dw.amp$p,
                                                      dw.res$r, dw.res$dw, dw.res$p,
                                                      rss, res1[,k][1], res1[,k][2])
  names(values) <- c(names(par$params), paste0(names(dw.amp)[1:3], "-amp"),
                     paste0(names(dw.res)[1:3], "-res"), c("rss", "ct", "eff"))
  rownames(values) <- c() #getting rid of arbitrary row names
  #plotting amplification curve
  xs = listdf$Cycle
if(plot){
  if(est$name == "l4"){
    plot(x=xs, y=l4_model(xs, b=par$params$b[k], c=par$params$c[k],
                          d=par$params$d[k], e=par$params$e[k]), type="l",  
         xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                                                   names(unlist(listdf)))])])), col = k, xaxt = "n")
    points(x=xs, y=listdf[,k+1], cex=0.45)
    legend("topleft", c(names(listdf)[k+1]), 
           col=k, lty=1, cex=0.65)
    #adding box around CT values (+/- 2 cycles)
    points(x=res1[,k][1], y=l4_model(res1[,k][1], b=par$params$b[k], 
                                     c=par$params$c[k], d=par$params$d[k], 
                                     e=par$params$e[k]), cex=0.8, pch=16) #CT point
    #boundaries for vertical lines
    clip(min(xs), max(xs), 
         l4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k]),
         l4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k])) 
    abline(v=res1[,k][1]-2, lty=1, col=k) ; abline(v=res1[,k][1]+2, lty=1, col=k) 
    #boundaries for horizontal lines
    clip(res1[,k][1]-2, res1[,k][1]+2, min(listdf[ ,2:length(listdf)][[k]]),
                                       max(listdf[ ,2:length(listdf)][[k]])) 
    abline(h=l4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                      d=par$params$d[k], e=par$params$e[k]), lty=1, col=k)
    abline(h=l4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                      d=par$params$d[k], e=par$params$e[k]), lty=1, col=k)
    #filling in the +/- squares for CT value
    polygon(x = c(res1[,k][1]-2, res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, res1[,k][1]-2), 
            y = c(l4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  l4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  l4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  l4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  l4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k])), 
                  col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
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
    #adding box around CT values (+/- 2 cycles)
    points(x=res1[,k][1], y=b4_model(res1[,k][1], b=par$params$b[k], 
                                     c=par$params$c[k], d=par$params$d[k], 
                                     e=par$params$e[k]), cex=0.8, pch=16) #CT point
    #boundaries for vertical lines
    clip(min(xs), max(xs), 
         b4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k]),
         b4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k])) 
    abline(v=res1[,k][1]-2, lty=1, col=k) ; abline(v=res1[,k][1]+2, lty=1, col=k) 
    #boundaries for horizontal lines
    clip(res1[,k][1]-2, res1[,k][1]+2, min(listdf[ ,2:length(listdf)][[k]]),
                                       max(listdf[ ,2:length(listdf)][[k]])) 
    abline(h=b4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                     d=par$params$d[k], e=par$params$e[k]), lty=1, col=k)
    abline(h=b4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                     d=par$params$d[k], e=par$params$e[k]), lty=1, col=k)
    #filling in the +/- squares for CT value
    polygon(x = c(res1[,k][1]-2, res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, res1[,k][1]-2), 
            y = c(b4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  b4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  b4_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  b4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  b4_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k])), 
            col= rgb(0,0,0,alpha=0.15))
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
    #adding box around CT values (+/- 2 cycles)
    points(x=res1[,k][1], y=l5_model(res1[,k][1], b=par$params$b[k], #CT point
                                                  c=par$params$c[k], d=par$params$d[k], 
                                                  e=par$params$e[k]), cex=0.8, pch=16) 
    #boundaries for vertical lines
    clip(min(xs), max(xs), 
         l5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k]),
         l5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k])) 
    abline(v=res1[,k][1]-2, lty=1, col=46) ; abline(v=res1[,k][1]+2, lty=1, col=k) 
    #boundaries for horizontal lines
    clip(res1[,k][1]-2, res1[,k][1]+2, min(listdf[ ,2:length(listdf)][[k]]),
                                       max(listdf[ ,2:length(listdf)][[k]])) 
    abline(h=l5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                     d=par$params$d[k], e=par$params$e[k]), lty=1, col=k)
    abline(h=l5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                     d=par$params$d[k], e=par$params$e[k]), lty=1, col=k)
    #filling in the +/- squares for CT value
    polygon(x = c(res1[,k][1]-2, res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, res1[,k][1]-2), 
            y = c(l5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  l5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  l5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  l5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  l5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k])), 
            col= rgb(0,0,0,alpha=0.15))
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
    
    #adding box around CT values (+/- 2 cycles)
    points(x=res1[,k][1], y=b5_model(res1[,k][1], b=par$params$b[k], #CT point
                                                  c=par$params$c[k], d=par$params$d[k], 
                                                  e=par$params$e[k]), cex=0.8, pch=16) 
    #boundaries for vertical lines
    clip(min(xs), max(xs), 
         b5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k]),
         b5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k])) 
    abline(v=res1[,k][1]-2, lty=1, col=46) ; abline(v=res1[,k][1]+2, lty=1, col=k) 
    #boundaries for horizontal lines
    clip(res1[,k][1]-2, res1[,k][1]+2, min(listdf[ ,2:length(listdf)][[k]]),
                                       max(listdf[ ,2:length(listdf)][[k]])) 
    abline(h=b5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                     d=par$params$d[k], e=par$params$e[k]), lty=1, col=k)
    abline(h=b5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                     d=par$params$d[k], e=par$params$e[k]), lty=1, col=k) 
    #filling in the +/- squares for CT value
    polygon(x = c(res1[,k][1]-2, res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, res1[,k][1]-2), 
            y = c(b5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  b5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  b5_model(res1[,k][1]+2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]), 
                  b5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k]),
                  b5_model(res1[,k][1]-2, b=par$params$b[k], c=par$params$c[k],
                                          d=par$params$d[k], e=par$params$e[k])), 
            col= rgb(0,0,0,alpha=0.15))
  }
  #plotting residuals
  plot(y = resids[[k]][1:length(listdf$Cycle)], 
       x = par$fits[[k]]$DATA$Cycles[1:length(listdf$Cycle)], 
       ylim = range(unlist(resids)), 
       xlab = "Cycle", ylab="Fluorescence Residual")
  abline(h=0) ; points(x=res1[,k][1], y=0, cex=0.8, pch=17)
  clip(min(xs), max(xs), 
       min(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]), 
       max(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]))
  abline(v=res1[,k][1]-2, lty=1, col=1) ; abline(v=res1[,k][1]+2, lty=1, col=1)
  clip(res1[,k][1]-2, res1[,k][1]+2, min(resids[[k]]), max(resids[[k]]))
  abline(h=min(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]), lty=1, col=1)
  abline(h=max(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]), lty=1, col=1)
  #filling in the +/- squares for CT value for resids
  polygon(x = c(res1[,k][1]-2, res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, res1[,k][1]-2), 
          y = c(min(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]), 
                max(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]),
                max(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]),
                min(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)]),
                min(resids[[k]][round(res1[,k][1]-2, digits=0):round(res1[,k][1]+2, digits=0)])),
          col= rgb(0,0,0,alpha=0.15))
}
  #return matrix of values
  return(values)
}
