plot_subset <- function(xs, params, params2, params3, params4){
  plot(x=xs, y=l4(xs, b=params$b[1], c=params$c[1],
                      d=params$d[1], e=params$e[1]), type="l",  
       xlab = "Cycle", ylab = "Fluorescence", ylim=c(min(df[2:ncol(df)]), max(df)))
  for(k in 1:nrow(params)){
    lines(x=xs, y=l4(xs, b=params$b[k], c=params$c[k],
                     d=params$d[k], e=params$e[k]), col=1)
    lines(x=xs, y=l5(xs, b=params2$b[k], c=params2$c[k],
                     d=params2$d[k], e=params2$e[k], f=params2$f[k]), col=2)
    lines(x=xs, y=b4(xs, b=params3$b[k], c=params3$c[k],
                     d=params3$d[k], e=params3$e[k]), col=3)
    lines(x=xs, y=b5(xs, b=params4$b[k], c=params4$c[k],
                     d=params4$d[k], e=params4$e[k], f=params4$f[k]), col=4)

    for(j in 2:ncol(A)){
      for(h in 1:nrow(params)){
        points(x=xs, y=listdf[[h]][,j], cex=0.45)
      }
    }
  }
  legend("topleft", c("b5", "b4", "l5", "l4"), 
         col=1:4, ncol=4, lty=1, cex=0.5)
}

plot_subset(Cycle, DF1, DF2, DF3, DF4)
