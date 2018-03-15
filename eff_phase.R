genparamsbase <- function(df, par){
subsetsb <- subsets
for (i in 1:8){
  for (k in 2:13){
    subsetsb[[i]][[k]] <- df[[i]][[k]] - par$params$c[[i]]
    }
}
  return(list(subsetsb))
}

plot_eff <- function(xs, listdf){
  for (k in 1:length(listdf)){
    for (i in 2:13){
        if (i == 2){
          plot(x = xs, y = listdf[[k]][[i]], type = "p", col = 1, xlab = "Cycle", ylab = "Fluorescence Efficiency")
        }
          points(x = xs, y = listdf[[k]][[i]], col = i-1)
        for (j in 1:39){
          points(x = xs[j], y = rowMeans(listdf[[k]][j,-1]), pch = 16)
        }
    }
    title(main = names(listdf[k]))
    
    }
  }