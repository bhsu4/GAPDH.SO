
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/")
df <- read.csv(file="GAPDH.SO.csv", header = TRUE, sep = ",")
Cycle = c(1:40)
library(data.table)
subsets <- lapply(LETTERS[1:8], function(k) cbind(Cycle,df[,colnames(df) %like% k]))
names(subsets) <- LETTERS[1:8]
#list2env(subsets, envir=.GlobalEnv)
#C <- C[,-c(1)]

library(qpcR)
source("GAPDH.SO/genparams.R")
df_b5 <- genparams(est=b5, listdf=subsets)

## to get residuals
resids <- lapply(df_b5$fits, resid)

source("GAPDH.SO/b5.R")
plot_b5(Cycle, subsets, df_b5$params)

for(k in 1:12){
  if(k == 1) plot(y=resids[[1]][1:40], x=df_b5$fits[[1]]$DATA$Cycles[1:40], ylim=range(resids[[1]]), type="l")
  if(k > 1){
    ind2 <- 40*k
    ind1 <- ind2-39
    lines(y=resids[[1]][ind1:ind2], x=df_b5$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
  }
}


