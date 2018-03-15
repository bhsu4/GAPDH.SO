genparamsbase <- function(df, par){
subsetsb <- subsets
subsetsbE <- subsets
for (i in 1:8){
  for (k in 2:13){
    subsetsb[[i]][[k]] <- df[[i]][[k]] - par$params$c[[i]]
    }
}
for (i in 1:8){
  for (k in 2:13){
    for (j in 2:40){
    subsetsbE[[i]][[k]][j-1,] <- (subsetsb[[i]][[k]][j,])/(subsetsb[[i]][[k]][j-1,])
    }
  }
}
  return(list(subsetsb, subsetsbE))
}

subsetsb_b5 <- genparamsbase(subsets, df_b5)


testing[[1]][[1]] <- (subsetsb_b5[[1]][[1]][2,]/subsetsb_b5[[1]][[1]][1,])
                     
#  []-params[,"c"][i])/(j[k-1,]-params[,"c"][i])

subsetsb[[1]][[1]] <- subsets [[1]][[2]] - df_b5$params$c[[1]]


for (h in 1:8){
assign(paste0("sub", h), 
       for (i in 2:8){
         for (k in 2:13){
         lapply(listdf, function(j) j[,k]-params[,"c"][i])
         }
       }
      )
}

lapply(listdf, function(j) j[,k]-params[,"c"][i])




lapply(listdf, function(j) 
  for (i in 2:8){
      for (k in 2:13){
          listEn[i] <- j[,k]-params[,"c"][i]
  }
})

#### test list df
lapply(listdf, function(j) 
  for (k in 2:8){
    for (i in 1:8){
      listEn[i] <- (j[k,]-params[,"c"][i])/(j[k-1,]-params[,"c"][i])
    }
  })   

      listEn[1] <- (A[2,]-params[,"c"][i])/(A[])
      
for (i in 1:8){

for(k in 2:13){
    if (k<3){
      for (i in 1:8){
        assign(paste0("subs", i))
          data.frame(c(listdf[[i]][,k]-DF1[,"c"][i]))
      }
          data.frame(c(listdf[[i]][,k]-DF1[,"c"][i]))
    }

subsparams(DF1)      
      

### figuring out assign

for (i in 1:2){
  
  assign(paste0("testingassign", i)
    
     data.frame(1,3,2)     
     data.frame(c(2, 4, 6))
         
         
         )
  }




### using previous
if(any(gsub("[[:alpha:]]","", result[[1]]$MODEL$name) == "4") == "TRUE") {
  for (k in 1:n){
    if (k < 2) {
      params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
      test <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                           params[,"e"]))
    }
    params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
    test[,k] <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                             params[,"e"]))
  }
  colnames(test) <- c("A", "B", "C", "D", "E", "F", "G", "H")
  newtest <- data.frame(t(test))
  h = which(listmod == result[[1]]$MODEL$name)
  assign(paste0("DF", h), newtest)
}
















listEn <- paste(LETTERS[1:8], "En", sep="")
result <- lapply(listdf, function(j) 
  for (k in 2:nrow(A)){
    for(i in 1:8)
      listEn[i] <- j[k,]/j[k-1,]
  })

m <- data.frame(matrix(0, ncol = 13, nrow = 1))

BEn <- B[k,]/B[k-1,]
CEn <- C[k,]/C[k-1,]
DEn <- D[k,]/D[k-1,]
EEn <- E[k,]/E[k-1,]
FEn <- F[k,]/F[k-1,]
GEn <- G[k,]/G[k-1,]
HEn <- H[k,]/H[k-1,]

for (k in 2:nrow(A)){
  AEn[k-1,] <- A[k,]/A[k-1,]
  BEn[k-1,] <- B[k,]/B[k-1,]
  CEn[k-1,] <- C[k,]/C[k-1,]
  DEn[k-1,] <- D[k,]/D[k-1,]
  EEn[k-1,] <- E[k,]/E[k-1,]
  FEn[k-1,] <- F[k,]/F[k-1,]
  GEn[k-1,] <- G[k,]/G[k-1,]
  HEn[k-1,] <- H[k,]/H[k-1,]
}

library(plyr)
colnames(m) <- colnames(A)
AEn <- rbind(m, AEn)
colnames(m) <- colnames(B)
BEn <- rbind(m, BEn)
colnames(m) <- colnames(C)
CEn <- rbind(m, CEn)
colnames(m) <- colnames(D)
DEn <- rbind(m, DEn)
colnames(m) <- colnames(E)
EEn <- rbind(m, EEn)
colnames(m) <- colnames(F)
FEn <- rbind(m, FEn)
colnames(m) <- colnames(G)
GEn <- rbind(m, GEn)
colnames(m) <- colnames(H)
HEn <- rbind(m, HEn)

Enlisted <- list(AEn, BEn, CEn, DEn, EEn, FEn, GEn, HEn)
Enlisted[[1]][,2:13]

CycleEn = 1:39
xcs <- CycleEn
plot_subset(CycleEn, Enlisted)

dev.off()
plot(x=1:39, y=Enlisted[[1]][,2], 
     xlab = "Cycle", ylab = "Efficiency", col=1)
for(j in (1:8)){
  if (j==1){
    for(k in 3:13){
      points(x=1:39, y=Enlisted[[j]][,k], col=1)
    }
  }
  for(k in 2:13) {
    points(x=1:39, y=Enlisted[[j]][,k], col=j)
  }
}

Enlisted2 <- list(AEn, BEn, CEn, DEn, EEn, FEn, GEn, HEn)
plot(x=listdf[[1]]$A10, y=Enlisted2[[1]][,2])




