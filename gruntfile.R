setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/")
df <- read.csv(file="GAPDH.SO.csv", header = TRUE, sep = ",")
Cycle = c(1:40)
library(data.table)
subsets <- lapply(LETTERS[1:8], function(k) cbind(Cycle,df[,colnames(df) %like% k]))
names(subsets) <- LETTERS[1:8]
#list2env(subsets, envir=.GlobalEnv)
subsets$C <- subsets$C[,-c(1)]

library(qpcR)
source("GAPDH.SO/genparams.R")
df_b5 <- genparams(est=b5, listdf=subsets)
df_l5 <- genparams(est=l5, listdf=subsets)
df_b4 <- genparams(est=b4, listdf=subsets)
df_l4 <- genparams(est=l4, listdf=subsets)

## to get residuals
source("GAPDH.SO/plot_resid.R")
b5resids <- lapply(df_b5$fits, resid)
plot_resid(df_b5, b5resids)

source("GAPDH.SO/b5_model.R")
plot_b5(Cycle, subsets, df_b5)
source("GAPDH.SO/b4_model.R")
plot_b4(Cycle, subsets, df_b4)
source("GAPDH.SO/l4_model.R")
plot_l4(Cycle, subsets, df_l4)
source("GAPDH.SO/l5_model.R")
plot_l5(Cycle, subsets, df_l5)

source("GAPDH.SO/plot_subset.R") #by model
plot_subset(Cycle, df_l4, df_l5, df_b4, df_b5, subsets)
source("GAPDH.SO/plot_together.R") 
plot_together(Cycle, df_l4, df_l5, df_b4, df_b5, subsets)

### subsetting with efficiencies / run fluo eff plots

source("GAPDH.SO/eff_phase.R") 
subsetsb_b5 <- genparamsbase(subsets, df_b5, subsets)

#subsetsbE <- lapply(1, matrix, data= NA, nrow=39, ncol=13)

#subsetsbE <- subsets
#for (i in 1:8){
#  for (k in 2:13){
#    for (j in 2:40){
#      subsetsbE[[i]][[k]] <- (subsetsb_b5[[i]][[k]][j,])/(subsetsb_b5[[i]][[k]][j-1,])
#    }
#  }
#}

#testing[[1]][[1]] <- (subsetsb_b5[[1]][[1]][2,]/subsetsb_b5[[1]][[1]][1,])

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
testingsubset <- lapply(LETTERS[1:8], function(k) cbind(Cycles, testing2[,colnames(testing2) %like% k]))
names(testingsubset) <- LETTERS[1:8]
testingsubset$C <- testingsubset$C[,-c(1)]

plot_efftest(Cycles, df_b5, testingsubset)


listdfb <- subsetb_b5
for(k in 1:length(listdfb)){ #turn it into 39 rows
  listdfb[[k]] <- subsetsb_b5[[k]][-40, ]
}
for(k in 1:length(listdfb)){
  for(j in 2:(nrow(listdfb$A)+1)){ #efficiency big, check baseline sub fluo is non-neg without division OK
    listdfb[[k]][j-1,] <- (subsetsb_b5[[k]][j,])/(subsetsb_b5[[k]][j-1,])
  }
}


source("GAPDH.SO/eff_phase.R")
Cycles = 1:39
subsetsb_b5 <- genparamsbase(subsets, df_b5, subsets)
plot_eff(Cycles, subsetsb_b5) #data's efficiency curve

curve_b5 <- curvefunc(df_b5, "b5_model", subsetsb_b5)
curveb_b5 <- genparamsbase(curve_b5, df_b5, subsets)
plot_eff(Cycles, curveb_b5)

# Plotting Double-Log Function

eff_b5 <- eff_values(subsetsb_b5)

leff_b5 <- log(log(eff_b5))

subsetscyc <- subsets
for(k in 1:length(subsetscyc)){ #turn it into 39 rows
  subsetscyc[[k]] <- subsets[[k]][-40, ]
}
plot(log(subsetscyc[[1]][,2]), log(eff_b5[[1]][,2]))
points(log(subsetscyc[[1]][,4]), log(eff_b5[[1]][,4]), col = 2)
points(log(subsetscyc[[1]][,5]), log(eff_b5[[1]][,5]), col = 3)
points(log(subsetscyc[[1]][,6]), log(eff_b5[[1]][,6]), col = 4)
points(log(subsetscyc[[1]][,7]), log(eff_b5[[1]][,7]), col = 5)

subsets_df <- data.frame(subsets)
cc = 1:11
hi <- lm(subsets_df[1:11,2] ~ cc, data = subsets_df)
hi$coefficients

hi <- lm(subsets[[1]][1:11,2] ~ cc, data = subsets)
listparams[[1]][1] <- list(hi$coefficients)

result <- lapply(subsets, function(k) lm(k ~ cc))

params222[[1]][,1] <- hi$coefficients

#creating ground phase slanted baseline model
params222 <- subsetsb_b5
for(k in 1:length(subsetsb_b5)){ #turn it into 39 rows
  params222[[k]] <- subsetsb_b5[[k]][-(3:40), -1]
}
cc <- 1:14
for (i in 1:8){
  for (k in 2:13){
    hi <- lm(subsets[[i]][1:14, k] ~ cc, data = subsets)
    params222[[i]][,k-1] <- hi$coefficients
    
  }
}


subsets_slb5 <- eff_b5

for (i in 1:8){
  for (k in 2:13){
    for (j in 2:40){
    subsets_slb5[[i]][[k]][[j-1]] <- (subsets[[i]][[k]][[j]]-(params222[[i]][k-1][1,]+params222[[i]][k-1][2,]*j))/(subsets[[i]][[k]][[j-1]]-(params222[[i]][k-1][1,]+params222[[i]][k-1][2,]*(j-1)))
    }
  }
}

plot_eff(Cycles, subsets_slb5)
