
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
subsetsb_b5 <- genparamsbase(subsets, df_b5)

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
subsetstest <- lapply(LETTERS[1:8], function(k) cbind(Cycles, testing2[,colnames(testing2) %like% k]))
names(subsetstest) <- LETTERS[1:8]
subsetstest$C <- subsetstest$C[,-c(1)]

plot_eff(Cycles, subsetstest)
