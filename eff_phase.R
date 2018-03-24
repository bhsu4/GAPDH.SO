genparamsbase <- function(df, par, sub){
subsetsb <- subsets
for (i in 1:8){
  for (k in 2:13){
    subsetsb[[i]][[k]] <- df[[i]][[k]] - par$params$c[[i]]
    }
}
  return(subsetsb)
}

curvefunc <- function(par, modelname, listdf){
  curve_values <- listdf
  for (k in 1:8){
    for (j in 2:13){
      for (i in 1:40){
        curve_values[[k]][[i,j]] <- do.call(what=modelname, args=list(x=i, b=par$params$b[[k]], c=par$params$c[[k]],
                                                                      d=par$params$d[[k]], e=par$params$e[[k]], f=par$params$f[[k]]))
      }
    }
  }
  return(curve_values)
}

plot_eff <- function(xs, listdf){
  
  listdfb <- listdf
  for(k in 1:length(listdfb)){ #turn it into 39 rows
    listdfb[[k]] <- listdf[[k]][-40, ]
  }
  for(k in 1:length(listdfb)){
    for(j in 2:(nrow(listdfb$A)+1)){ #efficiency big, check baseline sub fluo is non-neg without division OK
      listdfb[[k]][j-1,] <- (listdf[[k]][j,])/(listdf[[k]][j-1,])
    }
  }
  listdfb$C <- listdfb$C[,-c(1)]
  
  for (k in 1:length(listdfb)){
    for (i in 2:length(listdfb$A)){
      if (i == 2){
        plot(x = xs, y = listdfb[[k]][,i], type = "p", col = 1, xlab = "Cycle", ylab = "Fluorescence Efficiency",
             ylim=c(0,5))
      }
      points(x = xs, y = listdfb[[k]][,i], col = i-1)
      for (j in 1:nrow(listdfb$A)){
        points(x = xs[j], y = rowMeans(listdfb[[k]][j,-1]), pch = 19)
      }
    }
    title(main = names(listdfb[k]))
  }
}
