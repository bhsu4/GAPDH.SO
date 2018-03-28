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

plot_eff <- function(xs, listdf, sub, sl_df, baseline, slantk){
  
  #xs = Cycles
  #listdf = chosen model by genbaseparams (ex. subsetsb_b5)
  #sub = data used for regression for slant
  #sl_df = non-baselined df for slant method
  #baseline = flat or slant
  #slantk = cycle number for end of ground phase
  
  if(baseline == "flat") {
    
    listdfb <- listdf
    for(k in 1:length(listdfb)){ #turn it into 39 rows
      listdfb[[k]] <- listdf[[k]][-40, ]
    }
    for(k in 1:length(listdfb)){
      for(j in 2:(nrow(listdfb$A)+1)){ #efficiency big, check baseline sub fluo is non-neg without division OK
        listdfb[[k]][j-1,] <- (listdf[[k]][j,])/(listdf[[k]][j-1,])
      }
    }
    for (k in 1:length(listdfb)){
      for (i in 2:length(listdfb$A)){
        if (i == 2){
          plot(x = xs, y = listdfb[[k]][,i], type = "p", col = 1, xlab = "Cycle", ylab = "Fluorescence Efficiency",
               ylim = c(0,5))
        }
        points(x = xs, y = listdfb[[k]][,i], col = i-1)
        for (j in 1:nrow(listdfb$A)){
          points(x = xs[j], y = rowMeans(listdfb[[k]][j,-1]), pch = 19)
        }
      }
      title(main = names(listdfb[k]))
    }
  }
  
  if(baseline == "slant") {
    
    breg <- listdf
    for(k in 1:length(listdf)){ #turn it into 2 rows
      breg[[k]] <- listdf[[k]][-(3:(nrow(listdf$A))), -1]
    }
    cc <- 1:slantk
    for (i in 1:length(listdf)){
      for (k in 2:ncol(listdf$A)){
        regmod <- lm(sub[[i]][1:slantk, k] ~ cc, data = sub)
        breg[[i]][,k-1] <- regmod$coefficients
      }
    }
    subsets_sl <- listdf
    for(k in 1:length(listdf)){ #turn it into 39 rows
      subsets_sl[[k]] <- listdf[[k]][-40, ]
    }
    for (i in 1:8){
      for (k in 2:13){
        for (j in 2:40){
          subsets_sl[[i]][[k]][[j-1]] <- (sl_df[[i]][[k]][[j]]-(breg[[i]][k-1][1,]+ breg[[i]][k-1][2,]*j))/(sl_df[[i]][[k]][[j-1]]-(breg[[i]][k-1][1,]+breg[[i]][k-1][2,]*(j-1)))
        }
      }
    }
    for (k in 1:length(subsets_sl)){
      for (i in 2:length(subsets_sl$A)){
        if (i == 2){
          plot(x = xs, y = subsets_sl[[k]][,i], type = "p", col = 1, xlab = "Cycle", ylab = "Fluorescence Efficiency",
               ylim = c(0,5))
        }
        points(x = xs, y = subsets_sl[[k]][,i], col = i-1)
        for (j in 1:nrow(listdf$A)){
          points(x = xs[j], y = rowMeans(subsets_sl[[k]][j,-1]), pch = 19)
        }
      }
      title(main = names(subsets_sl[k]))
    }  
  }
}

eff_values <- function(listdf){
listdfb <- listdf
for(k in 1:length(listdfb)){ #turn it into 39 rows
  listdfb[[k]] <- listdf[[k]][-40, ]
}
for(k in 1:length(listdfb)){
  for(j in 2:(nrow(listdfb$A)+1)){ #efficiency big, check baseline sub fluo is non-neg without division OK
    listdfb[[k]][j-1,] <- (listdf[[k]][j,])/(listdf[[k]][j-1,])
  }
}
return(listdfb)
}