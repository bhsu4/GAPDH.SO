par(xpd = TRUE, mar = par()$mar + c(0,0,.5,0))

l5 <- function(x, b, c, d, e, f) {
  c + ( (d-c) / (1+exp(b*(log(x)-log(e))))^f )
}

plot_l5 <- function(xs, params){
  plot(x=xs, y=l5(xs, b=params$b[1], c=params$c[1],
                  d=params$d[1], e=params$e[1], f=params$f[1]), type="l",  
       xlab="Cycle", ylab="Fluorescence", ylim=c(min(df[2:ncol(df)]), max(df)))
  for(k in 2:nrow(params)){
    lines(x=xs, y=l5(xs, b=params$b[k], c=params$c[k],
                     d=params$d[k], e=params$e[k], f=params$f[k]),  
          col=k)
    for(j in 2:ncol(A)){
      for(h in 1:nrow(params)){
        points(x=xs, y=listdf[[h]][,j], cex=0.45)
      }
    }
  }
  legend(0, y=max(df)+100000, c(LETTERS[1:nrow(params)]), 
         col=1:nrow(params), ncol=8, lty=1, cex=0.5)
}

