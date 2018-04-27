b5_model <- function(x, b, c, d, e, f) {
  c+(d-c)/(1+exp(b*(x-e)))^f
}

plot_b5 <- function(xs, listdf, par){
  plot(x=xs, y=b5_model(xs, b=par$params$b[1], c=par$params$c[1],
                  d=par$params$d[1], e=par$params$e[1], f=par$params$f[1]), type="l",  
       xlab="Cycle", ylab="Fluorescence", ylim=c(min(df[2:ncol(df)]), max(df)))
  for(k in 2:length(subsets)){
    lines(x=xs, y=b5_model(xs, b=par$params$b[k], c=par$params$c[k],
                     d=par$params$d[k], e=par$params$e[k], f=par$params$f[k]),  
          col=k)
    for(j in 2:ncol(listdf$A)){
      for(h in 1:length(subsets)){
        points(x=xs, y=listdf[[h]][,j], cex=0.45)
      }
    }
  }
  legend("topleft", c(LETTERS[1:length(subsets)]), 
         col=1:length(subsets), ncol=2, lty=1, cex=0.5)
}

subplot_b5 <- function(xs, listdf, par){
  plot(x=xs, y=b5_model(xs, b=par$params$b[1], c=par$params$c[1],
                        d=par$params$d[1], e=par$params$e[1], f=par$params$f[1]), type="l",  
       xlab="Cycle", ylab="Fluorescence", ylim=c(min(df[2:ncol(df)]), max(df)))
  for(k in 2:length(listdf)){
    lines(x=xs, y=b5_model(xs, b=par$params$b[k], c=par$params$c[k],
                           d=par$params$d[k], e=par$params$e[k], f=par$params$f[k]),  
          col=k)
    for(j in 2:ncol(listdf)){
      points(x=xs, y=listdf[,j], cex=0.45)
    }
  }
  legend("topleft", c(names(listdf)[2:13]), 
         col=1:length(listdf), ncol=2, lty=1, cex=0.5)
}