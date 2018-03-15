genparamsbase <- function(df, par){
subsetsb <- subsets
for (i in 1:8){
  for (k in 2:13){
    subsetsb[[i]][[k]] <- df[[i]][[k]] - par$params$c[[i]]
    }
}
  return(list(subsetsb, subsetsbE))
}

subsetsb_b5 <- genparamsbase(subsets, df_b5)


subsetsbE <- lapply(1, matrix, data= NA, nrow=39, ncol=13)

subsetsbE <- subsets
for (i in 1:8){
  for (k in 2:13){
    for (j in 2:40){
      subsetsbE[[i]][[k]] <- (subsetsb_b5[[i]][[k]][j,])/(subsetsb_b5[[i]][[k]][j-1,])
    }
  }
}

testing[[1]][[1]] <- (subsetsb_b5[[1]][[1]][2,]/subsetsb_b5[[1]][[1]][1,])


testing <- data.frame(subsetsb_b5)
testing2 <- testing
for (k in 2:40){
testing2[k-1,] <- testing[k,]/testing[k-1,]
}
testing2 <- abs(testing2[-40,])
testing2 <- testing2[!grepl("Cycle",colnames(testing2))]
Cycles <- 1:39
testing2 <- cbind(Cycles, testing2)                       
colnames(testing2) <- colnames(df)
subsetstest <- lapply(LETTERS[1:8], function(k) cbind(Cycles, testing2[,colnames(testing2) %like% k]))
names(subsetstest) <- LETTERS[1:8]
subsetstest$C <- subsetstest$C[,-c(1)]


plot(x = Cycles, y = subsetstest[[1]][[2]], type="p", col = 1, ylab = "Fluorescence Efficiency")
for (i in 3:13){
  for (j in 1:39){
    
  points(x = Cycles, y = subsetstest[[1]][[i]], col = i-1)
  points(x = Cycles, y = rowMeans(subsetstest[[1]][j,-1]), pch  = 16 )
    
  }
}

for (j in 1:39){
ind2 <- 39*j
ind1 <- ind2-38
dftest[[i]][ind1] <- rowMeans(subsetstest[[1]][j,-1])
}


plot_eff <- function(xs, listdf){
for (k in 1:length(listdf)){
  for (i in 2:13){
    for (j in 1:39){
    if (i == 2){
    plot(x = xs, y = listdf[[k]][[i]], type = "p", col = 1, xlab = "Cycle", ylab = "Fluorescence Efficiency")
    }
    points(x = xs, y = listdf[[k]][[i]], col = i-1)
    points(x = xs[j], y = rowMeans(listdf[[k]][j,-1]), pch = 16)
  }
}
  title(main = names(listdf[k]))
}
}

plot_eff(Cycles, subsetstest)