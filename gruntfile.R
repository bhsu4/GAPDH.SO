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
Cycles = 1:39
subsetsb_b5 <- genparamsbase(subsets, df_b5, subsets)

#creating ground phase slanted baseline model
plot_eff(Cycles, subsetsb_b5, subsets, subsets, "slant", 11)
plot_eff(Cycles, subsetsb_b5, subsets, subsets, "flat", 11)
plot_eff(Cycles, subsetsb_b5, subsets, curve_b5, "slant", 11)

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

#### confirmation on curve values fitted

testconf <- function(x, b, e, f){
  ((1+exp(b*((x-1)-e)))/(1+exp(b*(x-e))))^f
}

plot_testconf <- function(xs, listdf, par){
  plot(x=xs, y=testconf(xs, b=par$params$b[1], e=par$params$e[1], f=par$params$f[1]), type="l",  
       xlab="Cycle", ylab="Fluorescence")
  for(k in 2:length(subsets)){
    lines(x=xs, y=testconf(xs, b=par$params$b[k], e=par$params$e[k], f=par$params$f[k]), col=k)
  }
  legend("topright", c(LETTERS[1:length(subsets)]), 
         col=1:length(subsets), ncol=2, lty=1, cex=0.5)
}

plot_testconf(Cycles, subsets, df_b5) #looks like eff. from curve values


