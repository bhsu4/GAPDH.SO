df <- read.csv(file="GAPDH.SO.csv", header = TRUE, sep = ",")
Cycle = c(1:40)
library(data.table)
subsets <- lapply(LETTERS[1:8], function(k) cbind(Cycle,df[,colnames(df) %like% k]))
names(subsets) <- LETTERS[1:8]
list2env(subsets, envir=.GlobalEnv)
C <- C[,-c(1)]