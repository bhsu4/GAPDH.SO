###Actual v. Slant Efficiency of GAPDH.SO

#xs = Cycles
#listdf = chosen model by genbaseparams (ex. subsetsb_b5)
#sub = data used for regression for slant
#sl_df = non-baselined df for slant method
#baseline = flat or slant
#slantk = cycle number for end of ground phase

#no axis labels
plot_eff_noaxis <- function(xs, listdf, sub, sl_df, baseline, slantk){

  if(baseline == "flat") {
    listdfb <- listdf[-40,]
    
    for(j in 2:40){ #efficiency big, check baseline sub fluo is non-neg without division OK
      listdfb[j-1,] <- (listdf[j,])/(listdf[j-1,])
    }
    for (i in 2:13){
      if (i == 2){
        plot(x = 1:39, y = listdfb[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5), yaxt="n", xaxt="n")
      }
      points(x = 1:39, y = listdfb[,i], col = i-1)
      for (j in 1:nrow(listdfb)){
        points(x = Cycles[[j]], y = rowMeans(listdfb[j,-1]), pch = 19, yaxt="n", xaxt = "n")
      }
    }
  }
  
  #listdf = subsetsb_b5[[1]]  
  if(baseline == "slant") {
    breg <- listdf[-(3:40), -1]
    cc <- 1:slantk
    for (k in 2:13){
      regmod <- lm(sub[1:11, k] ~ cc, data = subsets)
      breg[,k-1] <- regmod$coefficients
    }
    subsets_sl <- listdf[-40, ]
    #sl_df = subsets[[1]]  
    for (k in 2:13){
      for (j in 2:40){
        subsets_sl[[k]][[j-1]] <- (sl_df[[k]][[j]]-(breg[k-1][1,]+ breg[k-1][2,]*j))/(sl_df[[k]][[j-1]]-(breg[k-1][1,]+breg[k-1][2,]*(j-1)))
      }
    }
    
    for (i in 2:13){
      if (i == 2){
        plot(x = Cycles, y = subsets_sl[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5), yaxt="n", xaxt="n")
      }
      points(x = Cycles, y = subsets_sl[,i], col = i-1)
      for (j in 1:nrow(listdf)){
        points(x = Cycles[j], y = rowMeans(subsets_sl[j,-1]), pch = 19, yaxt = "n", xaxt="n")
      }
    }
  }  
}

#slant actual values y and x axis
plot_eff_yxact <- function(xs, listdf, sub, sl_df, baseline, slantk){
  
  if(baseline == "flat") {
    listdfb <- listdf[-40,]
    
    for(j in 2:40){ #efficiency big, check baseline sub fluo is non-neg without division OK
      listdfb[j-1,] <- (listdf[j,])/(listdf[j-1,])
    }
    for (i in 2:13){
      if (i == 2){
        plot(x = 1:39, y = listdfb[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5), yaxt="n", xaxt="n")
      }
      points(x = 1:39, y = listdfb[,i], col = i-1)
      for (j in 1:nrow(listdfb)){
        points(x = Cycles[[j]], y = rowMeans(listdfb[j,-1]), pch = 19, yaxt="n", xaxt = "n")
      }
    }
  }
  
  #listdf = subsetsb_b5[[1]]  
  if(baseline == "slant") {
    breg <- listdf[-(3:40), -1]
    cc <- 1:slantk
    for (k in 2:13){
      regmod <- lm(sub[1:11, k] ~ cc, data = subsets)
      breg[,k-1] <- regmod$coefficients
    }
    subsets_sl <- listdf[-40, ]
    #sl_df = subsets[[1]]  
    for (k in 2:13){
      for (j in 2:40){
        subsets_sl[[k]][[j-1]] <- (sl_df[[k]][[j]]-(breg[k-1][1,]+ breg[k-1][2,]*j))/(sl_df[[k]][[j-1]]-(breg[k-1][1,]+breg[k-1][2,]*(j-1)))
      }
    }
    
    for (i in 2:13){
      if (i == 2){
        plot(x = Cycles, y = subsets_sl[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5))
      }
      points(x = Cycles, y = subsets_sl[,i], col = i-1)
      for (j in 1:nrow(listdf)){
        points(x = Cycles[j], y = rowMeans(subsets_sl[j,-1]), pch = 19)
      }
    }
  }  
}

#slant fitted values onlyy x axis
plot_eff_xfit <- function(xs, listdf, sub, sl_df, baseline, slantk){
  
  if(baseline == "flat") {
    listdfb <- listdf[-40,]
    
    for(j in 2:40){ #efficiency big, check baseline sub fluo is non-neg without division OK
      listdfb[j-1,] <- (listdf[j,])/(listdf[j-1,])
    }
    for (i in 2:13){
      if (i == 2){
        plot(x = 1:39, y = listdfb[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5), yaxt="n", xaxt="n")
      }
      points(x = 1:39, y = listdfb[,i], col = i-1)
      for (j in 1:nrow(listdfb)){
        points(x = Cycles[[j]], y = rowMeans(listdfb[j,-1]), pch = 19, yaxt="n", xaxt = "n")
      }
    }
  }
  
  #listdf = subsetsb_b5[[1]]  
  if(baseline == "slant") {
    breg <- listdf[-(3:40), -1]
    cc <- 1:slantk
    for (k in 2:13){
      regmod <- lm(sub[1:11, k] ~ cc, data = subsets)
      breg[,k-1] <- regmod$coefficients
    }
    subsets_sl <- listdf[-40, ]
    #sl_df = subsets[[1]]  
    for (k in 2:13){
      for (j in 2:40){
        subsets_sl[[k]][[j-1]] <- (sl_df[[k]][[j]]-(breg[k-1][1,]+ breg[k-1][2,]*j))/(sl_df[[k]][[j-1]]-(breg[k-1][1,]+breg[k-1][2,]*(j-1)))
      }
    }
    
    for (i in 2:13){
      if (i == 2){
        plot(x = Cycles, y = subsets_sl[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5), yaxt = "n")
      }
      points(x = Cycles, y = subsets_sl[,i], col = i-1)
      for (j in 1:nrow(listdf)){
        points(x = Cycles[j], y = rowMeans(subsets_sl[j,-1]), pch = 19, yaxt = "n")
      }
    }
  }  
}

#slant fitted values onlyy x axis
plot_eff_yaxt <- function(xs, listdf, sub, sl_df, baseline, slantk){
  
  if(baseline == "flat") {
    listdfb <- listdf[-40,]
    
    for(j in 2:40){ #efficiency big, check baseline sub fluo is non-neg without division OK
      listdfb[j-1,] <- (listdf[j,])/(listdf[j-1,])
    }
    for (i in 2:13){
      if (i == 2){
        plot(x = 1:39, y = listdfb[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5), yaxt="n", xaxt="n")
      }
      points(x = 1:39, y = listdfb[,i], col = i-1)
      for (j in 1:nrow(listdfb)){
        points(x = Cycles[[j]], y = rowMeans(listdfb[j,-1]), pch = 19, yaxt="n", xaxt = "n")
      }
    }
  }
  
  #listdf = subsetsb_b5[[1]]  
  if(baseline == "slant") {
    breg <- listdf[-(3:40), -1]
    cc <- 1:slantk
    for (k in 2:13){
      regmod <- lm(sub[1:11, k] ~ cc, data = subsets)
      breg[,k-1] <- regmod$coefficients
    }
    subsets_sl <- listdf[-40, ]
    #sl_df = subsets[[1]]  
    for (k in 2:13){
      for (j in 2:40){
        subsets_sl[[k]][[j-1]] <- (sl_df[[k]][[j]]-(breg[k-1][1,]+ breg[k-1][2,]*j))/(sl_df[[k]][[j-1]]-(breg[k-1][1,]+breg[k-1][2,]*(j-1)))
      }
    }
    
    for (i in 2:13){
      if (i == 2){
        plot(x = Cycles, y = subsets_sl[,i], type = "p", col = 1, xlab = "Cycle", 
             ylab = "Fluorescence Efficiency", ylim = c(0,5), xaxt = "n")
      }
      points(x = Cycles, y = subsets_sl[,i], col = i-1)
      for (j in 1:nrow(listdf)){
        points(x = Cycles[j], y = rowMeans(subsets_sl[j,-1]), pch = 19, xaxt = "n")
      }
    }
  }  
}

dev.off()
par(oma=c(4,4,0.5,0.5),mar=c(0.25,0.25,0,0),mfrow=c(4,4),pch=16)

subsetsb_b4 <- genparamsbase(subsets, df_b4, subsets)
curve_b4 <- curvefunc(df_b4, b4_model, subsetsb_b4, 4)
plot_eff_noaxis(Cycles, subsetsb_b4[[1]], subsets[[1]], subsets[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_b4[[1]], subsets[[1]], subsets[[1]], "flat", 11)
plot_eff_noaxis(Cycles, subsetsb_b4[[1]], subsets[[1]], curve_b4[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_b4[[1]], subsets[[1]], curve_b4[[1]], "flat", 11)

subsetsb_b5 <- genparamsbase(subsets, df_b5, subsets)
curve_b5 <- curvefunc(df_b5, b5_model, subsetsb_b5, 5)
plot_eff_yaxt(Cycles, subsetsb_b5[[1]], subsets[[1]], subsets[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_b5[[1]], subsets[[1]], subsets[[1]], "flat", 11)
plot_eff_noaxis(Cycles, subsetsb_b5[[1]], subsets[[1]], curve_b5[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_b5[[1]], subsets[[1]], curve_b5[[1]], "flat", 11)

subsetsb_l4 <- genparamsbase(subsets, df_l4, subsets)
curve_l4 <- curvefunc(df_l4, l4_model, subsetsb_l4, 4)
plot_eff_noaxis(Cycles, subsetsb_l4[[1]], subsets[[1]], subsets[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_l4[[1]], subsets[[1]], subsets[[1]], "flat", 11)
plot_eff_noaxis(Cycles, subsetsb_l4[[1]], subsets[[1]], curve_l4[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_l4[[1]], subsets[[1]], curve_l4[[1]], "flat", 11)

subsetsb_l5 <- genparamsbase(subsets, df_l5, subsets)
curve_l5 <- curvefunc(df_l5, l5_model, subsetsb_l5, 5)
plot_eff_yxact(Cycles, subsetsb_l5[[1]], subsets[[1]], subsets[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_l5[[1]], subsets[[1]], subsets[[1]], "flat", 11)
plot_eff_xfit(Cycles, subsetsb_l5[[1]], subsets[[1]], curve_l5[[1]], "slant", 11)
plot_eff_noaxis(Cycles, subsetsb_l5[[1]], subsets[[1]], curve_l5[[1]], "flat", 11)

mtext(text="Cycles", side=1, line=2, outer=TRUE)
mtext(text="Efficiency", side=2, line=2, outer=TRUE)

##Residuals of GAPDH.SO with branching
## to get residuals
source("GAPDH.SO/plot_resid.R")
b5resids <- lapply(df_b5$fits, resid)

dev.off()
par(oma=c(2,2,0.5,0.5),mar=c(1.25,2,0.75,0.75),mfrow=c(2,4),pch=16)
plot_resid(subsets, df_b5)

plot_resid_axis <- function(params){
  
  for (i in 1:length(params$fits)){
    if(i == 1 ){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-40000,40000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
        axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    } 
    if(i == 2 |i==3 | i==4 ){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-40000,40000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
        #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    }
    if(i == 5){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-40000,40000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
        axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    }
    
    if(i==6 | i==7 | i==8){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-40000,40000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
        #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    }
    title(main= paste(names(params$fits[i]), sub=params$fits$A$MODEL$name, sep = ", "))
  }
}
b5resids <- lapply(df_b5$fits, resid)

dev.off()
m <- matrix(c(1,2,3,4,5,6,7,8),nrow = 2, ncol = 4, byrow = TRUE)
layout(mat = m)
par(oma=c(4,4,4,4),mar=c(0.5,0.75,1,0))

plot_resid_axis(df_b5)
mtext(text="Residual", side=2, line=2, outer=T)
mtext(text="Cycle", side=1, line=2, outer=T)

##AutoCorrelation of GAPDH.SO l5 Model 

#colors
my_colors <- brewer.pal(nlevels(as.factor(l5dat$ind)), "Set1")
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
my_colors <- add.alpha(my_colors, alpha=0.75) #adding transparency

#scatterplotMatrix (slow) 
l5dat[,"logrss"] <- log10(l5dat$rss)
hist(l5dat$logrss, n=30, xlim = c(2.5,5))
l5dat[,"ind"] <- ifelse(l5dat$logrss <= 3.3, "small", ifelse(l5dat$logrss >3.7, "big", "ok"))

spm(~dw.res+ct | ind, data=l5dat, oma=c(4,4,6,4), reg.line = F, smooth = F, cex=0.75, 
    pch=c(17,18,19,20), col=my_colors, main="Scatter plot l5", legend.plot=F)
par(xpd=TRUE)
legend("right", as.vector(unique(l5dat$ind)[1:3]), fill=my_colors)

#contour plot of scatterplot matrix
library(ggplot2)
ggplot(l5dat,  mapping = aes(dw.res, ct, col = ind, lty = ind)) +
  geom_density2d(contour = TRUE, size=0.75, linemitre = 3, bins=5) + xlim(range(l5dat$dw.res)) + ylim(range(l5dat$ct)) +
  scale_color_brewer(palette="Reds") + theme_bw()
#key point here!
#smaller dw.res (more ac) with larger RSS. there's something bad about 
#sig model that goes beyond the eq., increasing variance

#miRcompData Branching Residuals
try_b5 <- genparams(est=b5, listdf=try)

plot_resid_axis <- function(params){
  cyclength = lapply(params$fits, function(x) max(x$DATA[,1]))
  
  for (i in 1:length(params$fits)){
    if(i == 1 ){
      for(k in 1:4){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:cyclength[[i]]], 
                        x=params$fits[[i]]$DATA$Cycles[1:cyclength[[i]]], 
                        ylim=c(-200,200), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
        axis(2,at=seq(-200,200,50)) # add a new x-axis
        if(k > 1){
          ind2 <- cyclength[[i]]*k
          ind1 <- ind2-(cyclength[[i]]-1)
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[i]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    } 
    if(i == 2 |i==3 | i==4 | i==5){
      for(k in 1:4){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:cyclength[[i]]], 
                        x=params$fits[[i]]$DATA$Cycles[1:cyclength[[i]]], 
                        ylim=c(-200,200), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
        #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- cyclength[[i]]*k
          ind1 <- ind2-(cyclength[[i]]-1)
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[i]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    }
    if(i == 6){
      for(k in 1:4){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:cyclength[[i]]], 
                        x=params$fits[[i]]$DATA$Cycles[1:cyclength[[i]]], 
                        ylim=c(-200,200), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
        axis(2,at=seq(-200,200,50)) # add a new x-axis
        if(k > 1){
          ind2 <- cyclength[[i]]*k
          ind1 <- ind2-(cyclength[[i]]-1)
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[i]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    }
    
    if(i==7 | i==8 | i==9 | i==10){
      for(k in 1:4){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:cyclength[[i]]], 
                        x=params$fits[[i]]$DATA$Cycles[1:cyclength[[i]]], 
                        ylim=c(-200,200), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
        #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- cyclength[[i]]*k
          ind1 <- ind2-(cyclength[[i]]-1)
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[i]]$DATA$Cycles[ind1:ind2], col=k)
        }
      }
    }
    title(main= paste(names(params$fits[i]), sub=params$fits$A$MODEL$name, sep = ", "))
  }
}
b5resids <- lapply(try_b5$fits, resid)
dev.off()
m <- matrix(c(1,2,3,4,5,6,7,8,9,10),nrow = 2, ncol = 5, byrow = TRUE)
layout(mat = m)
par(oma=c(4,4,4,4),mar=c(0.5,0.75,1,0))

plot_b5(try, try_b5) #plot of the data fluo
subplot_b5(try)
plot_resid_axis(try_b5) #each rep vs. 1 uniform curve
#for fits that have mediocre, good signal fits, we see the same branching trend
#fits with noise don't have that pattern

#residuals with each replicate as curve of its own
b5_resultsG <- sub_genparams(b5, subsets$G)
b4_resultsG <- sub_genparams(b4, subsets$G)
l5_resultsG <- sub_genparams(l5, subsets$G)
l4_resultsG <- sub_genparams(l4, subsets$G)

subplot_resid_4x <- function(params){
  
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1) plot(y=resids[[k]][1:40], 
                    x=params$fits[[1]]$DATA$Cycles[1:40], 
                    type="l", ylim = c(-11900,9000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
                    xlab="", ylab="", yaxt = "n")
    #axis(side=2, at=seq(-2000, 2000, by=2000)) #-8000,8000 and -2000,2000
    
    if(k > 1){
      ind2 <- 40*k
      ind1 <- ind2-39*(k)-(k-1)
      lines(y=resids[[k]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
}

subplot_resid_4xy <- function(params){
  
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1) plot(y=resids[[k]][1:40], 
                    x=params$fits[[1]]$DATA$Cycles[1:40], 
                    type="l", ylim = c(-11900,9000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
                    xlab="", ylab="", yaxt = "n")
    axis(side=2, at=seq(-8000, 8000, by=4000)) #-8000,8000 and -2000,2000
    
    if(k > 1){
      ind2 <- 40*k
      ind1 <- ind2-39*(k)-(k-1)
      lines(y=resids[[k]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
}

subplot_resid_5xy <- function(params){
  
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1) plot(y=resids[[k]][1:40], 
                    x=params$fits[[1]]$DATA$Cycles[1:40], 
                    type="l", ylim = c(-3000,3000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
                    xlab="", ylab="", yaxt = "n", xaxt="n")
    axis(side=2, at=seq(-2000, 2000, by=2000)) #-8000,8000 and -2000,2000
    
    if(k > 1){
      ind1 = 1; ind2=40
      #ind2 <- 40*k
      #ind1 <- ind2-39*(k)-(k-1)
      lines(y=resids[[k]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
}

subplot_resid_5x <- function(params){
  
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1) plot(y=resids[[k]][1:40], 
                    x=params$fits[[1]]$DATA$Cycles[1:40], 
                    type="l", ylim = c(-3000,3000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
                    xlab="", ylab="", yaxt = "n", xaxt="n")
    #axis(side=2, at=seq(-2000, 2000, by=2000)) #-8000,8000 and -2000,2000
    
    if(k > 1){
      ind1 = 1 ; ind2 = 40
      #ind2 <- 40*k
      #ind1 <- ind2-39*(k)-(k-1)
      lines(y=resids[[k]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
}

dev.off()
par(mfrow=c(2,2))
par(oma=c(4,4,0.5,4),mar=c(0,0.25,0,0))

subplot_resid_5xy(l5_resultsG)
subplot_resid_5x(b5_resultsG)
subplot_resid_4xy(l4_resultsG)
subplot_resid_4x(b4_resultsG)
mtext(text="Cycles", side=1, line=2, outer=TRUE)
mtext(text= "Residuals", side=2, line=2, outer=TRUE)

#miRcomp of each rep as own curve residual
subplot_resid_each <- function(tester){
  
b5_try <- sub_genparams(b5, tester)
b4_try <- sub_genparams(b4, tester)
l5_try <- sub_genparams(l5, tester)
l4_try <- sub_genparams(l4, tester)
cyclength = tester[,1]

tryresid <- function(x) {tryCatch({resid(x)}, error = function(e) {rep(NA, max(cyclength))})}
residsl5 <- lapply(l5_try$fits, tryresid)
residsb5 <- lapply(b5_try$fits, tryresid)
residsl4 <- lapply(l4_try$fits, tryresid)
residsb4 <- lapply(b4_try$fits, tryresid)

#range5 <- round_any(range(residsl5, residsb5, na.rm=TRUE), 100, f=ceiling)
#range4 <- round_any(range(residsl4, residsb4, na.rm=TRUE), 100, f=ceiling)
#if(any(abs(range(residsl5, residsb5, na.rm=TRUE)) < 100)==TRUE || any(abs(range(residsl4, residsb4, na.rm=TRUE)) < 100)==TRUE ){
#  range5 <- round_any(range(residsl5, residsb5, na.rm=TRUE), 25, f=ceiling)
#  range4 <- round_any(range(residsl4, residsb4, na.rm=TRUE), 25, f=ceiling)
#}

par(mfrow=c(2,2))
par(oma=c(4,4,0.5,4),mar=c(0.25,0.25,0,0))

for(x in list(residsl5, residsb5, residsl4, residsb4)){
  if(identical(x, residsl5) == TRUE){
  for(k in 1:4){
  if(k==1){
    plot(x=cyclength, y=x[[1]], col = k, type = "l", ylim = c(-150,150), 
         xlab="", ylab="", xaxt = "n", yaxt = "n")
    axis(side=2, at=seq(-100,100, by=50)) 
  }
    lines(x=cyclength, y=x[[k]], col = k)
    }
  }
  else if(identical(x, residsb5) == TRUE){
    for(k in 1:4){
      if(k==1){
        plot(x=cyclength, y=x[[1]], col = k, type = "l", ylim = c(-150,150), 
             xlab="", ylab="", xaxt = "n", yaxt = "n")
      }
      lines(x=cyclength, y=x[[k]], col = k)
    }
  }
  else if(identical(x, residsl4) == TRUE){
    for(k in 1:4){
      if(k==1){
        plot(x=cyclength, y=x[[1]], col = k, type = "l", ylim = c(-150,150), 
             xlab="", ylab="", yaxt = "n")
        axis(side=2, at=seq(-100,100, by=50)) 
      }
      lines(x=cyclength, y=x[[k]], col = k)
    }
  } 
  else if(identical(x, residsb4) == TRUE){
    for(k in 1:4){
      if(k==1){
        plot(x=cyclength, y=x[[1]], col = k, type = "l", ylim = c(-150,150), 
             xlab="", ylab="", yaxt = "n")
      }
      lines(x=cyclength, y=x[[k]], col = k)
    }
  }
}
mtext(text="Cycles", side=1, line=2, outer=TRUE)
mtext(text= "Residuals", side=2, line=2, outer=TRUE)
}
subplot_resid_each(try$J)


unlist.genparams <- function(test){
  tst.list <- list() ; listdf.tst <- list() ; repnames <- list()
  for(i in 1:length(test)){
    for(j in 1:length(test[[i]])){ #length = rep of secondary list
      mincyc <- min(unlist(lapply(test[[i]], function(k) max(k$Cycle)))) #min cyc (if diff)
      tst.list[[i]] <- sapply(test[[i]], function(x) x$Rn[1:mincyc]) #list of diff sampleIDs, but df secondary
      listdf.tst[[i]] <- as.data.frame(cbind(seq(1:mincyc), tst.list[[i]])) 
      repnames <- lapply(LETTERS[1:length(test)], paste0, 1:length(test[[i]])) #df colnames
      names(listdf.tst[[i]]) <- c("Cycle", repnames[[i]][1:length(test[[i]])])
    }
  }
  names(listdf.tst) <- LETTERS[1:length(test)]
  return(listdf.tst)
}

#using highest Rn instead of dRn to see if pattern same
miRcompData2.1 <- miRcompData2[order(miRcompData2$Rn, decreasing = TRUE),] #KW5_2
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets")
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-576-3p_002351.Rda")
try.big <- unlist.genparams(tst) #list of lists organized here
subplot_resid_each(try.big$E)

load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-500_002428.Rda")
try.big2 <- unlist.genparams(tst) #KW8
subplot_resid_each(try.big2$H)



#different signal fits
library(dynlm) ; library(car)


##good, bad, non-signals, unidentified signals (failed exp)

#nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams(tst)
source("GAPDH.SO/plot_sig.R")
plot_sig(l5, try.good, plot=TRUE) #plots shown

#not nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst)

#run this to see all curves
plot_sig(l5, try, plot=TRUE) #plots shown
plot_sig(l5, try.good, plot=TRUE) #plots shown
plot_sig(l4, try.good, plot=TRUE) #plots shown

##nice data set
#good fit (L5) J4, H4!
#good fit (L4) D1
#poor fit (L4) I4!, I3, I2, C2!

##not nice data set!
#failed (L5) exp I2!, J2, C3!
#nonsignals (L5) D2, D3, D4, C4!

par(oma=c(2,2,0.5,0.5),mar=c(1.25,2,0.75,0.75),mfrow=c(4,2),pch=16)
plot_sig(l5, try.good, macro=4, z=8, plot=TRUE) #goodfit nicedata l5
plot_sig(l4, try.good, macro=4, z=9, plot=TRUE) #poorfit nicedata l4
plot_sig(l5, try, macro=4, z=3, plot=TRUE) #nonsignal notnicedata l5
plot_sig(l5, try, macro=2, z=9, plot=TRUE) #failed notnicedata l5

#turn line color same
#turn off grey rectangle polgyon for >46 cycles 



