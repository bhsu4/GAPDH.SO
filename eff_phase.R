genparamsbase <- function(df, par, sub){
subsetsb <- subsets
for (i in 1:8){
  for (k in 2:13){
    subsetsb[[i]][[k]] <- df[[i]][[k]] - par$params$c[[i]]
    }
}
  return(subsetsb)
}

plot_efftest <- function(xs, par, listdf, listdfb){
  
  listdfb <- listdf
  for(k in 1:length(listdfb)){ #turn it into 39 rows
    listdfb[[k]] <- listdf[[k]][-40, ]
  }
  for(k in 1:length(listdfb)){
    for(j in 2:(nrow(listdfb$A)+1)){ #efficiency big, check baseline sub fluo is non-neg without division OK
      listdfb[[k]][j-1,] <- (listdf[[k]][j,]-par$params$c[[k]])/(listdf[[k]][j-1,]-par$params$c[[k]])
    }
  }
  
  for (k in 1:length(listdfb)){
    for (i in 2:length(listdfb$A)){
      if (i == 2){
        plot(x = xs, y = listdfb[[k]][,i], type = "p", col = 1, xlab = "Cycle", ylab = "Fluorescence Efficiency")
      }
      points(x = xs, y = listdfb[[k]][,i], col = i-1)
      for (j in 1:nrow(listdfb$A)){
        points(x = xs[j], y = rowMeans(listdfb[[k]][j,-1]), pch = 3)
      }
    }
    title(main = names(listdfb[k]))
  }
}
