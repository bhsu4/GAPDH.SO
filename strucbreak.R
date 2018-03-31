strbreak_chow <- function(xs, ys){
  struc_result <- ys
  break_result <- ys
  for(k in 1:length(ys)){ #turn it into 31 rows
    struc_result[[k]] <- ys[[k]][-c(1:9),-1]
    break_result[[k]] <- ys[[k]][-c(1:9),-1]
  }
  for(i in 1:length(ys)){ #1:8
    for(j in 2:length(ys$A)){ #2:13
      for(k in 10:nrow(ys$A)){ #10:40
        result <- sctest(ys[[i]][,j] ~ xs, type = "Chow", point = k-5)
        struc_result[[i]][,j-1][k-9] <- result$p.value
        dum.xs=rep(1,nrow(ys$A))
        dum.xs[x>=k-5]=0
        result2 <- summary(lm(ys[[i]][,j] ~ xs*dum.xs))
        
        break_result[[i]][,j-1][k-9] <- pf(result2$fstatistic[1], result2$fstatistic[2], 
                                           result2$fstatistic[3], lower.tail=F)
      }
    }
  }
  return(list(struc_result, break_result))
}

