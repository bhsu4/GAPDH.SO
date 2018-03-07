b5 <- function(x, b, c, d, e, f) {
  c+(d-c)/(1+exp(b*(x-e)))^f
}

plot_b5 <- function(xs, listdf, params){
  plot(x=xs, y=b5(xs, b=params$b[1], c=params$c[1],
                  d=params$d[1], e=params$e[1], f=params$f[1]), type="l",  
       xlab="Cycle", ylab="Fluorescence", ylim=c(min(df[2:ncol(df)]), max(df)))
  for(k in 2:nrow(params)){
    lines(x=xs, y=b5(xs, b=params$b[k], c=params$c[k],
                     d=params$d[k], e=params$e[k], f=params$f[k]),  
          col=k)
    for(j in 2:ncol(listdf$A)){
      for(h in 1:nrow(params)){
        points(x=xs, y=listdf[[h]][,j], cex=0.45)
      }
    }
  }
  legend("topleft", c(LETTERS[1:nrow(params)]), 
         col=1:nrow(params), ncol=2, lty=1, cex=0.5)
}

