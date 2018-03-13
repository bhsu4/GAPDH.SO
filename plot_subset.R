plot_subset <- function(xs, par, par2, par3, par4, listdf){
  plot(x=xs, y=l4_model(xs, b=par$params$b[1], c=par$params$c[1],
                      d=par$params$d[1], e=par$params$e[1]), type="l",  
       xlab = "Cycle", ylab = "Fluorescence", ylim=c(min(df[2:ncol(df)]), max(df)))
  for(k in 1:12){
    lines(x=xs, y=l4_model(xs, b=par$params$b[k], c=par$params$c[k],
                     d=par$params$d[k], e=par$params$e[k]), col=1)
    lines(x=xs, y=l5_model(xs, b=par2$params$b[k], c=par$params$c[k],
                     d=par2$params$d[k], e=par2$params$e[k], f=par2$params$f[k]), col=2)
    lines(x=xs, y=b4_model(xs, b=par3$params$b[k], c=par$params$c[k],
                     d=par3$params$d[k], e=par3$params$e[k]), col=3)
    lines(x=xs, y=b5_model(xs, b=par4$params$b[k], c=par$params$c[k],
                     d=par4$params$d[k], e=par4$params$e[k], f=par4$params$f[k]), col=4)
    for(j in 2:12){
      for(h in 1:8){
        points(x=xs, y=listdf[[h]][,j], cex=0.45)
      }
    }
  }
  legend("topleft", c("l4", "l5", "b4", "b5"), 
         col=1:4, ncol=4, lty=1, cex=0.5)
}