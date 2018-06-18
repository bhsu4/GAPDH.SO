l5_model <- function(x, b, c, d, e, f) {
  c + ( (d-c) / (1+exp(b*(log(x)-log(e))))^f )
}

plot_l5 <- function(listdf, par){
  cycunl = lapply(listdf, function(x) unlist(x$Cycle))
  cycmax <- max(unlist(cycunl)) ; xs = 1:cycmax
  plot(x=xs, y=l5_model(xs, b=par$params$b[1], c=par$params$c[1],
                            d=par$params$d[1], e=par$params$e[1], f=par$params$f[1]), 
       type="l", xlab="Cycle", ylab="Fluorescence", 
       ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                                                 names(unlist(listdf)))])])))
  for(k in 2:length(listdf)){
    lines(x=xs, y=l5_model(xs, b=par$params$b[k], c=par$params$c[k],
                               d=par$params$d[k], e=par$params$e[k], 
                               f=par$params$f[k]), col=k)
    for(j in 2:length(try[[i]])){
      for(h in 1:length(listdf)){
        points(x=xs, y=listdf[[h]][,j][1:cycmax], cex=0.45)
      }
    }
  }
  legend("topleft", c(LETTERS[1:length(listdf)]), 
         col=1:length(listdf), ncol=2, lty=1, cex=0.5)
}

subplot_l5 <- function(listdf){
  #get fitted models of each replicate within list
  par <- lapply(listdf, function(x) sub_genparams(l5, x)) #KW test results
  #plot the model
  for(i in 1:length(listdf)){
    xs = listdf$Cycle
    plot(x=xs, y=l5_model(xs, b=par$params$b[1], c=par$params$c[1],
                              d=par$params$d[1], e=par$params$e[1], 
                              f=par$params$f[1]), type="l",  
         xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                                                   names(unlist(listdf)))])]))) #xaxt="n", yaxt="n"
    for(k in 2:length(par$params)){
      lines(x=xs, y=l5_model(xs, b=par$params$b[k], c=par$params$c[k],
                                 d=par$params$d[k], e=par$params$e[k], 
                                 f=par$params$f[k]), col=k)
      for(j in 2:length(listdf)){
        points(x=xs, y=listdf[,j], cex=0.45, col = j-1)
      }
    }
    legend("topleft", c(names(listdf)[2:length(listdf)]), 
           col=1:(length(listdf)-1), ncol=2, lty=1, cex=0.65)
    #x.intersp=0.25, text.width=c(rep(0,6), rep(0.5,6))
  }
}