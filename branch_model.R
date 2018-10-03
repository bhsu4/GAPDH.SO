setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/")
df <- read.csv(file="GAPDH.SO.csv", header = TRUE, sep = ",")
Cycle = c(1:40)
library(data.table)
subsets <- lapply(LETTERS[1:8], function(k) cbind(Cycle,df[,colnames(df) %like% k]))
#colnames(subsets$F) <- c("Cycle", paste0("F", 1:12))
setnames(df, "F6.1", "F7") #changing one column name
names(subsets) <- LETTERS[1:8]
#list2env(subsets, envir=.GlobalEnv)
subsets$C <- subsets$C[,-c(1)]

library(qpcR)
source("GAPDH.SO/genparams.R")
df_b5 <- genparams(est=b5, listdf=subsets)
df_l5 <- genparams(est=l5, listdf=subsets)
df_b4 <- genparams(est=b4, listdf=subsets)
df_l4 <- genparams(est=l4, listdf=subsets)

source("GAPDH.SO/b4_model.R")
plot_b4(subsets, df_b4)
subplot_b4(subsets)

res1 <- list()
for(i in 1:8){
res1[[i]] <- getPar(df_b4$fits[[i]], type = "curve", cp = "cpD2", eff = "sliwin")
res1[[i]] <- res1[[i]][-2,]
}
b4CT <- ceiling(do.call(cbind, res1))
#start at 19, end at 26?

#Proposition 1
man.tilda <- (c/(r*d1))*sigma

sigma = pkn/mkn * ykn


nexp = 26
#Weighted cconditional least squares estimator ???
pkn.reps <- matrix(rep(NA, 4*8), nrow = 8)
ykn.reps <- matrix(rep(NA, 4*8), nrow = 8)
pkn.each <- rep(NA, 4) ; ykn.each <- rep(NA, 4)
for(i in 1:8){
  for(j in 2:5){
    texp = b4CT[[i]]
    pkn.num = sum(subsets[[i]][,j][1:nexp]) - subsets[[i]][,j][texp] - sum(subsets[[i]][,j][1:(nexp-1)])
    pkn.den = sum(subsets[[i]][,j][texp:nexp])
    #empty frames
    pkn.each[j-1] <- pkn.num/pkn.den
    ykn.each[j-1] <- sum(subsets[[i]][,j])
  }
  #fill in matrix
  pkn.reps[i,1:4] <- pkn.each
  ykn.reps[i,1:4] <- ykn.each
}
mkn.reps <- (1+pkn.reps)^(4+1)
mkn <- rowSums(mkn.reps)
ykn <- rowSums(ykn.reps)
pkn <- rowSums(pkn.reps)

sigma = pkn/mkn*ykn






  
  
  