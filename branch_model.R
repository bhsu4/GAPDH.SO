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
pkn.reps <- matrix(rep(NA, 12*8), nrow = 8)
ykn.reps <- matrix(rep(NA, 12*8), nrow = 8)
pkn.each <- rep(NA, 12) ; ykn.each <- rep(NA, 12)
for(i in 1:8){
  for(j in 2:13){
    texp = b4CT[[i]]
    pkn.num = sum(subsets[[i]][,j][1:nexp]) - subsets[[i]][,j][texp] - sum(subsets[[i]][,j][1:(nexp-1)])
    pkn.den = sum(subsets[[i]][,j][texp:nexp])
    #empty frames
    pkn.each[j-1] <- pkn.num/pkn.den
    ykn.each[j-1] <- sum(subsets[[i]][,j])
  }
  #fill in matrix
  pkn.reps[i,1:12] <- pkn.each
  ykn.reps[i,1:12] <- ykn.each
}
mkn.reps <- (1+pkn.reps)^(12+1)
mkn <- rowSums(mkn.reps)
ykn <- rowSums(ykn.reps)
pkn <- rowSums(pkn.reps)

sigma = pkn/mkn*ykn

abs_quant <- function(df_mod, subsets, dil = "reps"){
 ##key lengths defined
  replength = length(subsets[[1]])-1 
  sublength = length(subsets)
 ##finding start and end of exp phase 
  res1 <- list()
  for(i in 1:sublength){
    res1[[i]] <- getPar(df_mod$fits[[i]], type = "curve", cp = "cpD2", eff = "sliwin")
    res1[[i]] <- res1[[i]][-2,]
  }
  CTval <- ceiling(do.call(cbind, res1))
###Proposition #1: m_{an}
 ##empty matrices for replication results
  pkn.reps <- matrix(rep(NA, replength*sublength), nrow = sublength)
  ykn.reps <- matrix(rep(NA, replength*sublength), nrow = sublength)
  pkn.each <- rep(NA, replength) ; ykn.each <- rep(NA, replength)
 ##finding replication values
  #pkn = estimator of efficiency for kth reaction in exp phase
  #mkn = 1 + pkn
  #ykn = fluorescence sum in exp phase
  for(i in 1:sublength){
    for(j in 2:(replength+1)){
      texp = CTval[[i]]
      pkn.num = sum(subsets[[i]][,j][1:nexp]) - subsets[[i]][,j][texp] - sum(subsets[[i]][,j][1:(nexp-1)])
      pkn.den = sum(subsets[[i]][,j][texp:nexp])
      #empty frames
      pkn.each[j-1] <- pkn.num/pkn.den
      ykn.each[j-1] <- sum(subsets[[i]][,j])
    }
    #fill in matrix
    pkn.reps[i,1:replength] <- pkn.each
    ykn.reps[i,1:replength] <- ykn.each
  }
 ##accumulated amounts inside Sigma
  mkn.reps <- (1+pkn.reps)^(replength+1)
  mkn <- rowSums(mkn.reps)
  ykn <- rowSums(ykn.reps)
  pkn <- rowSums(pkn.reps)
  sigma = pkn*ykn/mkn
 ##constants outside Sigma
  if(dil = "reps"){
    denom = 1/replength
  }
  else{
    d1n = sum(dil) #dil is a vector b/c D1(n) = 1/r(n) * Sigma(k to r(n)) of dk
    denom = (1/(replength*dk1))
  }
  #return(list(m = mkn, y = ykn, p = pkn))
}

abs_quant(df_b4, subsets = subsets)


rel_quant <- function(df_mod, subsetsC, subsetsT){
  ##key lengths defined
  replength = length(subsets[[1]])-1 
  sublength = length(subsets)
  ##finding start and end of exp phase 
  res1 <- list()
  for(i in 1:sublength){
    res1[[i]] <- getPar(df_mod$fits[[i]], type = "curve", cp = "cpD2", eff = "sliwin")
    res1[[i]] <- res1[[i]][-2,]
  }
  CTval <- ceiling(do.call(cbind, res1))
  ###Proposition #1: m_{an}
  ##empty matrices for replication results
  pkn.repsC <- matrix(rep(NA, replength*sublength), nrow = sublength)
  pkn.repsT <- matrix(rep(NA, replength*sublength), nrow = sublength)
  ykn.repsC <- matrix(rep(NA, replength*sublength), nrow = sublength)
  ykn.repsT <- matrix(rep(NA, replength*sublength), nrow = sublength)
  pkn.eachC <- rep(NA, replength) ; pkn.eachT <- rep(NA, replength)
  ykn.eachC <- rep(NA, replength) ; ykn.eachT <- rep(NA, replength)
  ##finding replication values
  #pkn = estimator of efficiency for kth reaction in exp phase
  #mkn = 1 + pkn
  #ykn = fluorescence sum in exp phase
  for(i in 1:sublength){
    for(j in 2:(replength+1)){
      texp = CTval[[i]]
      pkn.numC = sum(subsetsC[[i]][,j][1:nexp]) - subsetsC[[i]][,j][texp] - sum(subsetsC[[i]][,j][1:(nexp-1)])
      pkn.denC = sum(subsetsC[[i]][,j][texp:nexp])
      pkn.numT = sum(subsetsT[[i]][,j][1:nexp]) - subsetsT[[i]][,j][texp] - sum(subsetsT[[i]][,j][1:(nexp-1)])
      pkn.denT = sum(subsetsT[[i]][,j][texp:nexp])
      #empty frames
      pkn.eachC[j-1] <- pkn.numC/pkn.denC ; pkn.eachT[j-1] <- pkn.numT/pkn.denT
      ykn.eachC[j-1] <- sum(subsetsC[[i]][,j]); ykn.eachT[j-1] <- sum(subsetsT[[i]][,j])
    }
    #fill in matrix
    pkn.repsC[i,1:replength] <- pkn.eachC ; pkn.repsT[i,1:replength] <- pkn.eachT
    ykn.repsC[i,1:replength] <- ykn.eachC ; ykn.repsT[i,1:replength] <- ykn.eachT
  }
  ##accumulated amounts inside Sigma
  mkn.repsC <- (1+pkn.repsC)^(replength+1) ; mkn.repsT <- (1+pkn.repsT)^(replength+1)
  mknC <- rowSums(mkn.repsC) ; mknT <- rowSums(mkn.repsT)
  yknC <- rowSums(ykn.repsC) ; yknT <- rowSums(ykn.repsT)
  pknC <- rowSums(pkn.repsC) ; pknT <- rowSums(pkn.repsT)
  sigmaC = pknC*yknC/mknC ; sigmaT = pknT*yknT/mknT
  ##constants outside Sigma
  manC = (1/replength)*pknC*yknC/mknC
  manT = (1/replength)*pknT*yknT/mknT
  #return(list(m = mkn, y = ykn, p = pkn))
  return(list(mC = manC, mT = manT))
}
  
  