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

dfb5.ml1 <- list() ; dfb5.res1 <- list()
for(i in 1:8){
  dfb5.ml1[[i]] <- modlist(df_b5$fits[[i]], model=l5)
  dfb5.res1[[i]] <- getPar(dfb5.ml1[[i]], type = "curve", cp = "cpD2", eff = "sliwin")
}

plot_resid_axis <- function(params){
  
  for (i in 1:length(params$fits)){
    if(i == 1 ){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-30000,30000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
        axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
          abline(v=dfb5.res1[[i]][1,], lty='dashed')
        }
      }
    } 
    if(i == 2 |i==3 | i==4 ){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-30000,30000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
        #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
          abline(v=dfb5.res1[[i]][1,], lty='dashed')
          
        }
      }
    }
    if(i == 5){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-30000,30000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
        axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
          abline(v=dfb5.res1[[i]][1,], lty='dashed')
          
        }
      }
    }
    
    if(i==6 | i==7 | i==8){
      for(k in 1:12){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:40], 
                        x=params$fits[[1]]$DATA$Cycles[1:40], 
                        ylim=c(-30000,30000), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
        #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
        if(k > 1){
          ind2 <- 40*k
          ind1 <- ind2-39
          lines(y=resids[[i]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
          abline(v=dfb5.res1[[i]][1,], lty='dashed')
          
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

#baseline subtracted gapdhso
subsets_base <- vector('list', 8) ; subsets2 <- list()
#baseline subtract GAPDH
for(i in 1:8){
  for(j in 2:13){
    subsets_base[[i]][[1]] <- subsets[[i]][,1]
    subsets_base[[i]][[j]] <- subsets[[i]][,j] - min(subsets[[i]][,j])
  }
  subsets2[[i]] <- data.frame(matrix(unlist(subsets_base[[i]]), ncol=13))
}
names(subsets2) <- LETTERS[1:8]
df_b5_base <- genparams(est=b5, listdf=subsets2)
plot_resid(subsets2, df_b5_base)

plot_resid_axis(df_b5_base)
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
  geom_density2d(contour = TRUE, size=0.75, linemitre = 3, bins=5) + #xlim(c(0,0.1)) +
  xlim(range(l5dat$dw.res)) + 
  ylim(range(l5dat$ct)) +
  scale_color_brewer(palette="Reds") + theme_bw()
#key point here!
#smaller dw.res (more ac) with larger RSS. there's something bad about 
#sig model that goes beyond the eq., increasing variance

ggplot(l5dat,  mapping = aes(dw.amp, ct, col = ind, lty = ind)) +
  geom_density2d(contour = TRUE, size=0.75, linemitre = 3, bins=5) + xlim(c(0.02,0.075)) +
  #xlim(range(l5dat$dw.amp)) + 
  ylim(range(l5dat$ct)) +
  scale_color_brewer(palette="Set1") + theme_bw()

smoothScatter(x=l5dat$dw.amp, y=l5dat$ct)

x1 = l5dat[l5dat$ind == 'small',]
x2 = l5dat[l5dat$ind == 'ok',]
x3 = l5dat[l5dat$ind == 'big',]

alpharamp<-function(c1,c2, alpha=128) {stopifnot(alpha>=0 & alpha<=256);function(n) paste(colorRampPalette(c(c1,c2))(n), format(as.hexmode(alpha), upper.case=T), sep="")}
smoothScatter(x=x1$dw.amp,y=x1$ct,nrpoints=length(x1),cex=3, colramp=alpharamp("white",blues9))
par(new=T)
smoothScatter(x=x2$dw.amp,y=x2$ct,nrpoints=length(x2),cex=3,colramp= alpharamp("white","red"), axes=F, ann=F)
par(new=T)
smoothScatter(x=x3$dw.amp,y=x3$ct,nrpoints=length(x2),cex=3,colramp= alpharamp("white","green"), axes=F, ann=F)



###contour plot with points
library(MASS)  # in case it is not already loaded 
set.seed(101)
n <- 1000
X <- mvrnorm(n, mu=c(.5,2.5), Sigma=matrix(c(1,.6,.6,1), ncol=2))

## some pretty colors
library(RColorBrewer)
k <- 12
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
l5dat.na <- l5dat[which(  (!is.na(l5dat$dw.amp)) & (!is.na(l5dat$ct))),]
z <- kde2d(x=l5dat.na$dw.amp, y=l5dat.na$ct, n=5000)

plot(x=l5dat.na$dw.amp, y=l5dat.na$ct, xlab="X label", ylab="Y label", pch=19, cex=.4, xlim = c(0.02, 0.25))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2.5)

abline(h=mean(X[,2]), v=mean(X[,1]), lwd=2)
legend("topleft", paste("R=", round(cor(X)[1,2],2)), bty="n")



#heat map in contour
par(oma=c(2,2,0.5,0.5),mar=c(1.25,2,0.75,0.75),mfrow=c(2,2),pch=16)
commonTheme = list(labs(color="Density",fill="Density",
                        x="Durbin-Watson Test Statistic",
                        y="CT Value (SDM)"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))
#dynamic linear model, dw.amp
plot1 <- ggplot(data=l5dat.na,aes(dw.amp,ct)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  #geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  xlim(c(0.025, 0.065))+
  guides(alpha="none") +
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme

#dw.comp
l5dat2.na <- l5dat2[which(  (!is.na(l5dat2$dw.comp)) & (!is.na(l5dat2$ct))),]

plot2 <- ggplot(data=l5dat2.na,aes(dw.comp, ct)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
 # geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  xlim(c(0,2.5))+ #ylim(c(5,40)) +
  guides(alpha="none") +
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme

#side by side
require(gridExtra)
grid.arrange(plot1, plot2, ncol=2)


#ljung-box
commonTheme = list(labs(color="Density",fill="Density",
                        x="Durbin-Watson Test Statistic",
                        y="Ljung-Box P-value"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))
plot3 <- ggplot(data=l5dat2.na,aes(dw.comp,boxlj.p)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  ylim(0,0.05)+
  xlim(c(0,2.5))+
  guides(alpha="none") + geom_hline(yintercept=0.05, lty='dashed')+
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme
commonTheme = list(labs(color="Density",fill="Density",
                        x="Ljung-Box P-value",
                        y="Pearson Correlation"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))
plot4 <- ggplot(data=l5dat2.na,aes(boxlj.p,pcor)) + 
  #stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  #scale_fill_continuous(low="green",high="red") +
  #geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  ylim(-1,1)+
  xlim(c(0,1))+
  guides(alpha="none") + geom_vline(xintercept=0.05, lty='dashed')+
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme

grid.arrange(plot3, plot4, ncol=2)



  #ljung-box
ggplot(data=l5dat2.na,aes(boxlj.p, ct)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  #geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  #ylim(0,0.5)+
  #xlim(c(0,0.0001))+
  guides(alpha="none") + geom_hline(yintercept=0.05, lty='dashed')+
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme

ggplot(data=l5dat2.na,aes(dw.comp,pcor)) + 
  #stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  #scale_fill_continuous(low="green",high="red") +
  #geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  ylim(-1,1)+
  xlim(c(0,4))+
  guides(alpha="none") +
  geom_point(shape=16, cex=0.05, alpha=0.25) + commonTheme


#miRcompData Branching Residuals
plot_resid_axis <- function(params){
  cyclength = lapply(params$fits, function(x) max(x$DATA[,1]))
  
  for (i in 1:length(params$fits)){
    if(i == 1 ){
      for(k in 1:4){
        resids <- lapply(params$fits, resid)
        
        if(k == 1) plot(y=resids[[i]][1:cyclength[[i]]], 
                        x=params$fits[[i]]$DATA$Cycles[1:cyclength[[i]]], 
                        ylim=c(-400,400), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
        axis(2,at=seq(-400,400,100)) # add a new x-axis
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
                        ylim=c(-400,400), type="l", 
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
                        ylim=c(-400,400), type="l", 
                        xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
        axis(2,at=seq(-400,400,100)) # add a new x-axis
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
                        ylim=c(-400,400), type="l", 
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
#dRn try
try_b5 <- genparams(est=b5, listdf=try)

#Rn nice dataset
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams(tst)
try.good <- unlist.genparams.Rn(tst)
try.good_b5 <- genparams(est=b5, listdf=try.good)
plot_resid_axis(try.good_b5) #each rep vs. 1 uniform curve
mtext("Cycle", side = 1, line= 2, outer = TRUE)
mtext("Residual", side = 2, line= 2, outer = TRUE)

dev.off()
for(i in 1:10){
  for(j in 2:5){
    if(j==2){
      plot(x=try.good[[i]][,1], y=try.good[[i]][,j], type = 'l', 
           col= j-1, ylim = range(unlist(try.good[[i]][,2:5])))
    }
      lines(x=try.good[[i]][,1], y=try.good[[i]][,j], col = j-1)
  }
}



plot_b5(try, try_b5) #plot of the data fluo
subplot_b5(try)
plot_resid_axis(try_b5) #each rep vs. 1 uniform curve



#using highest Rn instead of dRn to see if pattern same
miRcompData2.1 <- miRcompData2[order(miRcompData2$Rn, decreasing = TRUE),] #KW5_2
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets")
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-576-3p_002351.Rda")
try.big <- unlist.genparams(tst) #list of lists organized here
try.big <- unlist.genparams.Rn(tst)
try.big_b5 <- genparams(est=b5, listdf=try.big)
plot_resid_axis(try.big_b5)

load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-500_002428.Rda")
try.big2 <- unlist.genparams(tst) #KW8
try.big2 <- unlist.genparams.Rn(tst) #KW8
try.big2_b5 <- genparams(est=b5, listdf=try.big2)
plot_resid_axis(try.big2_b5)


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
                    xlab="", ylab="", yaxt = "n", xaxt = 'n')
    #axis(side=2, at=seq(-2000, 2000, by=2000)) #-8000,8000 and -2000,2000
    #axis(side=1, at=seq(10, 40, by=10))
    
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

#lag1
lstarmodnolog1 <- list()
for (j in 1:8) lstarmodnolog1[[j]] <- lapply(subsets[[j]][,2:13], function(f) try(tsDyn::lstar(f, m=1, d=1)))
#lag2
lstarmodnolog2 <- list()
for (j in 1:8) lstarmodnolog2[[j]] <- lapply(subsets[[j]][,2:13], function(f) try(tsDyn::lstar(f, m=2, d=1)))
lstarmodnolog2 <- lapply(lstarmodnolog2, function(x) lapply(x, function(x) c(rep(NA,2), x$residuals)))

dev.off()
m <- matrix(c(1,2,5,3,4,6),nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, widths=c(100, 100, 100))
par(oma=c(4,4,0.5,4),mar=c(0.25,0.25,0,0))

subplot_resid_5xy(l5_resultsG)
subplot_resid_5x(b5_resultsG)
subplot_resid_4xy(l4_resultsG)
subplot_resid_4x(b4_resultsG)

for(i in 1:12){
  if(i==1){
    plot(x=2:40, y=lstarmodnolog1[[7]][[i]]$residuals, type = 'l', 
         col = 1, ylim = c(range(lapply(lstarmodnolog1[[7]], function(x) unlist(x$residuals)))),
         xaxt = "n", yaxt = "n")
    axis(side=4, at=seq(-6000, 6000, by=3000))
    #axis(side=1, at=seq(0, 40, by=10))
    
  }
  else{
    lines(x=2:40, y=lstarmodnolog1[[7]][[i]]$residuals, col=i)
  }
}
  
for(i in 1:12){
  if(i==1){
    plot(x=1:40, y=lstarmodnolog2[[7]][[i]], type = 'l', 
         col = 1, ylim = c(range(lstarmodnolog2[[7]], na.rm=TRUE)), 
         xaxt = "n", yaxt = "n")
    axis(side=4, at=seq(-3000, 3000, by=1500))
    axis(side=1, at=seq(0, 40, by=10))
    
  }
  else{
    lines(x=1:40, y=lstarmodnolog2[[7]][[i]], col=i)
  }
}
mtext(text="Cycle", side=1, line=2, outer=TRUE)
mtext(text= "Residual", side=2, line=2, outer=TRUE)





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


unlist.genparams.Rn <- function(test){
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
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets")
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
#poor fit (L5) I3
#failed (L5) exp I2!, J2, C3!
#nonsignals (L5) D2, D3, D4, C4!

par(oma=c(2,2,0.5,0.5),mar=c(1.25,2,0.75,0.75),mfrow=c(4,2),pch=16)
plot_sig(l5, try.good, macro=4, z=8, plot=TRUE) #goodfit nicedata l5
plot_sig(l4, try.good, macro=4, z=9, plot=TRUE) #poorfit nicedata l4
plot_sig(l5, try, macro=4, z=3, plot=TRUE) #nonsignal notnicedata l5
plot_sig(l5, try, macro=2, z=9, plot=TRUE) #failed notnicedata l5

######changing plot_sig to exclude legend, and xaxis label

plot_sigmacro <- function(est, listdf, macro=0, z=NULL, plot=FALSE){
  
  #start of microscoping: specific w term evaluation
  
  #plotting amplification curve
  if(macro > 0){
    w = macro
    #re-state specified z-list
    listdf <- listdf[[z]]
    #getting parameter estimates
    par <- sub_genparams(est, listdf)
    #finding the residuals
    try_resid <- function(x) tryCatch({resid(x)},
                                      error = function(e) rep(NA, max(listdf$Cycle)))
    resids <- lapply(par$fits, try_resid)
    #finding the RSS
    rss <- sum(resids[[w]]^2)
    #finding the DW-stat
    reg.amp <- dynlm(listdf[,w+1] ~ listdf$Cycle)
    reg.res <- try(dynlm(resids[[w]] ~ listdf$Cycle), silent = TRUE)
    #DW-statistics for amp and resid
    dw.amp <- durbinWatsonTest(reg.amp) 
    dw.res <- tryCatch({
      durbinWatsonTest(reg.res)
    }, error = function(e) {
      return(list(r=NA, dw=NA, p=NA))
    })  #cannot run with NA's
    #finding the CT value
    ml1 <- modlist(listdf, model = est)
    res1 <- getPar(ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    #combining all into df
    values <- data.frame(apply(par$params[w,], c(1,2), as.numeric))
    values[,(length(values)+1):(length(values)+9)] <- c(dw.amp$r, dw.amp$dw, dw.amp$p,
                                                        dw.res$r, dw.res$dw, dw.res$p,
                                                        rss, res1[,w][1], res1[,w][2])
    names(values) <- c(names(par$params), paste0(names(dw.amp)[1:3], "-amp"),
                       paste0(names(dw.res)[1:3], "-res"), c("rss", "ct", "eff"))
    rownames(values) <- c() #getting rid of arbitrary row names  
    if(plot){
      xs = listdf$Cycle
      if(est$name == "l4"){ #if NAs, model cannot run
        try(plot(x=xs, y=l4_model(xs, b=par$params$b[w], c=par$params$c[w],
                                  d=par$params$d[w], e=par$params$e[w]), type="l",  
                 xlab="Cycle", ylab="Fluorescence", ylim = c(range(unlist(listdf[,w+1]))),
                 # ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                 #                              names(unlist(listdf)))])])), 
                 col = 1, xaxt = "n")) #col=w
        points(x=xs, y=listdf[,w+1], cex=0.45)
        #adding box around CT values (+/- 2 cycles)
        try(points(x=res1[,w][1], y=l4_model(res1[,w][1], b=par$params$b[w], 
                                             c=par$params$c[w], d=par$params$d[w], 
                                             e=par$params$e[w]), cex=0.8, pch=16)) #CT point
        #boundaries for vertical lines
        if( (is.na(res1[,w][[1]]) == "TRUE") || (res1[,w][[1]] <=2) || (res1[,w][[1]] > (max(listdf$Cycle) - 2))){
          print(paste0(LETTERS[z], w, " " , "no ct")) 
        } #can't draw box for those with <=2 CT or NA
        else{
          #filling in the +/- squares for CT value
          polygon(x = c(res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2, res1[,w][1]-2), 
                  y = c(min(par("usr")), min(par("usr")), 
                        max(par("usr")), max(par("usr")), min(par("usr"))),
                  col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
        }
      }
      if(est$name == "b4"){
        try(plot(x=xs, y=b4_model(xs, b=par$params$b[w], c=par$params$c[w],
                                  d=par$params$d[w], e=par$params$e[w]), type="l",  
                 xlab="Cycle", ylab="Fluorescence", 
                 #ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                 #                             names(unlist(listdf)))])])), 
                 ylim=c(range(unlist(listdf[,w+1]))),
                 col = 1, xaxt = "n")) #col=w
        points(x=xs, y=listdf[,w+1], cex=0.45)
        #adding box around CT values (+/- 2 cycles)
        try(points(x=res1[,w][1], y=b4_model(res1[,w][1], b=par$params$b[w], 
                                             c=par$params$c[w], d=par$params$d[w], 
                                             e=par$params$e[w]), cex=0.8, pch=16)) #CT point
        #filling in the +/- squares for CT value
        if(is.na(res1[,w][[1]]) == "FALSE"){
          polygon(x = c(res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2, res1[,w][1]-2), 
                  y = c(min(par("usr")), min(par("usr")), 
                        max(par("usr")), max(par("usr")), min(par("usr"))),
                  col= rgb(0,0,0,alpha=0.15))
        }
        else{
          print(paste0(LETTERS[z], w, " ", "no ct"))
        }
      }
      if(est$name == "l5"){
        try(plot(x=xs, y=l5_model(xs, b=par$params$b[w], c=par$params$c[w],
                                  d=par$params$d[w], e=par$params$e[w],
                                  f=par$params$f[w]), type="l",  
                 xlab="Cycle", ylab="Fluorescence", 
                 #ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                 #                             names(unlist(listdf)))])])), 
                 ylim=c(range(unlist(listdf[,w+1]))),
                 col = 1, xaxt = "n")) #col=w
        points(x=xs, y=listdf[,w+1], cex=0.45)
        #adding box around CT values (+/- 2 cycles)
        try(points(x=res1[,w][1], y=l5_model(res1[,w][1], b=par$params$b[w], #CT point
                                             c=par$params$c[w], d=par$params$d[w], 
                                             e=par$params$e[w], f=par$params$f[w]), cex=0.8, pch=16)) 
        #boundaries for vertical lines
        if( (is.na(res1[,w][[1]]) == "TRUE") || (res1[,w][[1]] <= 2) || (res1[,w][[1]] > (max(listdf$Cycle) - 2))){
          print(paste0(LETTERS[z], w, " " , "no ct"))
        }
        else{
          
          #filling in the +/- squares for CT value
          polygon(x = c(res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2, res1[,w][1]-2), 
                  y = c(min(par("usr")), min(par("usr")), 
                        max(par("usr")), max(par("usr")), min(par("usr"))),
                  col= rgb(0,0,0,alpha=0.15))
        }
      }
      if(est$name == "b5"){ 
        try(plot(x=xs, y=b5_model(xs, b=par$params$b[w], c=par$params$c[w],
                                  d=par$params$d[w], e=par$params$e[w],
                                  f=par$params$f[w]), type="l",  
                 xlab="Cycle", ylab="Fluorescence", ylim=c(range(unlist(listdf[,w+1]))),
                 #ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                 #                             names(unlist(listdf)))])])), 
                 col = 1, xaxt = "n")) #col=w
        points(x=xs, y=listdf[,w+1], cex=0.45)
        
        #adding box around CT values (+/- 2 cycles)
        try(points(x=res1[,w][1], y=b5_model(res1[,w][1], b=par$params$b[w], #CT point
                                             c=par$params$c[w], d=par$params$d[w], 
                                             e=par$params$e[w], f=par$params$f[w]), cex=0.8, pch=16)) 
        #boundaries for vertical lines
        if(is.na(res1[,w][[1]]) == "FALSE"){
          #filling in the +/- squares for CT value
          polygon(x = c(res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2, res1[,w][1]-2), 
                  y = c(min(par("usr")), min(par("usr")), 
                        max(par("usr")), max(par("usr")), min(par("usr"))),
                  col= rgb(0,0,0,alpha=0.15))
        }
        else{
          print(paste0(LETTERS[z], w, " " , "no ct"))
        }
      }
      #plotting residuals
      if(unique(is.na(resids[[w]])) == "FALSE"){
        plot(y = resids[[w]][1:length(listdf$Cycle)], 
             x = par$fits[[w]]$DATA$Cycles[1:length(listdf$Cycle)], 
             ylim = range(unlist(resids)[which(!is.na(unlist(resids)))]), 
             xlab = "Cycle", xaxt = "n", yaxt = "n", ylab="Fluorescence Residual")
        abline(h=0) ; points(x=res1[,w][1], y=0, cex=0.8, pch=17)
        #yaxis on right
        axis(4,at=seq(-40,20,10)) # add a new x-axis
        
      }
      else{
        plot(1, type="n", xlab="n", ylab="", xlim=c(0, sum(is.na(resids[[w]]))), ylim=c(0,1))
      }
      #boundaires for horizontal lines
      if( (res1[,w][[1]] <= 2) || (is.na(res1[,w][[1]])) == "TRUE" ||
          (res1[,w][[1]] > 38 & max(listdf$Cycle) == 40) ||
          (res1[,w][[1]] > 44 & max(listdf$Cycle) == 46) ){ 
        #recall: listdf is with primary list called already
        print("check")
      }
      else{
        polygon(x = c(res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2, res1[,w][1]-2), 
                y = c(min(par("usr")), min(par("usr")), 
                      max(par("usr")), max(par("usr")), min(par("usr"))), col= rgb(0,0,0,alpha=0.15))
      }
      
      #axis(1,at=seq(0,40,10)) # add a new x-axis
      #mtext(text="Residual", side=2, line=2, outer=T)
      #mtext(text="Cycle", side=1, line=2, outer=T)
      
    } #plot bracket
    return(values)
  } #macro bracket
}

par(oma=c(2,5,0,3), mar=c(0,0,1,2), mfrow=c(5,2),pch=16)
plot_sigmacro(l5, try.good, macro=4, z=8, plot=TRUE) #goodfit nicedata l5
plot_sigmacro(l4, try.good, macro=4, z=9, plot=TRUE) #poorfit nicedata l4
plot_sigmacro(l5, try.good, macro=4, z=9, plot=TRUE) #poorfit nicedata l5
plot_sigmacro(l5, try, macro=4, z=3, plot=TRUE) #nonsignal notnicedata l5
plot_sigmacro(l5, try, macro=2, z=9, plot=TRUE) #failed notnicedata l5
mtext(text="Cycle", side=1, line=1, outer=T)
mtext(text="Fluorescence", side=2, line=3, outer=T)
mtext(text="Residual", side=4, line=1, outer=T)


###Autocorrelation boxplots
#nicedataset
source("GAPDH.SO/gen_results.R")
modparams <- lapply(try.good, function(x) sub_genparams(listdf=x, est=l5))
try.sig <- gen_results(try, modparams)
lstarres <- vector("list", 10) #empty list of subs: ABCD,..etc
cyclength <- unlist(lapply(try.good, nrow))
for(i in 1:10){ #running LSTAR model
  for(j in 2:5){
    #results for LSTAR model 
    lstarres[[i]][[j-1]] <- tryCatch({
      tsDyn::lstar(try.good[[i]][,j], m=mdim, d=klag)}, #d = lag found through AIC
      error=function(e) list(fitted.values=rep(NA, (cyclength[[i]]-(klag*mdim))), 
                             residuals=rep(NA, (cyclength[[i]]-(klag*mdim))),
                             model.specific=list(coefficients = rep(NA, (4+2*mdim)))))
    #if error, output NA, (no fit) w/ fitted values and residual rep length times minus klag
  }
}
lstarres.fits <- lapply(lstarres, function(x) list(fits=x))
try.lstar <- gen_results(try.good, lstarres.fits)

res.tog <- data.frame(1:40, try.sig$amp.dw, try.lstar$amp.dw)
colnames(res.tog) <- c("Cycle", "Sig", "Lstar")
mdata <- melt(res.tog, id=c("Cycle"))
mdata[,"value"] <- round(as.numeric(mdata[,"value"]),2)

require(ggplot2)
ggplot(data = mdata, aes(x=Cycle, y=value)) + geom_boxplot(aes(fill=variable))

###AIC plt GAPDH.SO supplemental

for(i in 1:8){
  for(j in 1:13){
    ts_subsets[[i]][,j] <- ts(subsets[[i]][,j])
  }
}
#names(ts_subsets$F)[8] <- "F7" #replace F.6.1

source("GAPDH.SO/lag_gen.R")
library(dynlm)
library(Hmisc)
ts_subsets <- ts(subsets)
fftest <- lapply(LETTERS[1:8], paste0, 1:12) #list of subset names
fftry <- lapply(fftest, function(x) lapply(x, get_lag_formulae, n=15)) #dynlm formula in list
ff_fitstest <- list()
for(j in 1:8) ff_fitstest[[j]] <- lapply(fftry[[j]], function(x) lapply(x, 
                                                          function(f) dynlm(formula = as.formula(f), 
                                                              data=ts_subsets[[j]]))) #output results of dynlm
ff_fitstestaic <- list()
for(j in 1:8) ff_fitstestaic[[j]] <- sapply(ff_fitstest[[j]], 
                                            function(x) sapply(x,AIC)) #output results of dynlm

for(j in 1:8){
  for(i in 1:12){
    if(i==1 & j==1){
      plot(x=1:15, y=ff_fitstestaic[[j]][,i], type = "p", ylim = c(range(ff_fitstestaic)), 
           ylab = "AIC", xlab = "Lag Term")
      #axis(4)
      #mtext('AIC',4,line=2)
      #mtext('Lag Term', 1, line=2)
    }
    points(x=1:15, y=ff_fitstestaic[[j]][,i], col = j)
  } #plotting all replication sets' AIC
}
for(j in 1:8){
  lines(spline(x = 1:15, y=rowMeans(ff_fitstestaic[[j]])), pch = 19, cex = 3, col = j)
} #smoothing line



###plotting ljung-box and pearcor for GAPDH, and nice try.goodd
#gapdh.so ljung-box
gap.each <- list() ; 
for(i in 1:8){
gap.each[[i]] <- sub_genparams(l5, subsets[[i]])
}
gap.each.resids <- lapply(gap.each, function(x) lapply(x$fits, resid))

gap.each.box <- vector('list', 8) ; gap.each.dw <- vector('list', 8) ; peacor <- vector('list', 8)
for(i in 1:8){
  for(j in 1:12){
  gap.each.box[[i]][[j]] <- Box.test(gap.each.resids[[i]][[j]], lag=1, type='Ljung-Box')
  gap.each.dw[[i]][[j]] <- sum((gap.each.resids[[i]][[j]]-Lag(gap.each.resids[[i]][[j]], 1))^2, na.rm=TRUE)/sum(gap.each.resids[[i]][[j]]^2)
  peacor[[i]][[j]] <- cor(gap.each.resids[[i]][[j]][2:40], 
                          gap.each.resids[[i]][[j]][1:39])
  }
}
gapdh.ljdw <- data.frame(matrix(c(unlist(gap.each.dw), unlist(lapply(gap.each.box, function(x) lapply(x, function(x) unlist(x$p.value)))), unlist(peacor)), ncol=3))
colnames(gapdh.ljdw) <- c("dw", "lj.p", "peacor")

par(mfrow=c(1,1))
plot(x=gapdh.ljdw$lj.p, y=gapdh.ljdw$peacor)

gapdh.ljdw.melt <- melt(gapdh.ljdw, id.vars=c())

par(mfrow=c(1,3))
gplot1 <- ggplot(data = gapdh.ljdw, aes(x="", y=dw)) + geom_boxplot() +xlab("") +ylab("")
gplot2<- ggplot(data = gapdh.ljdw, aes(x="", y=lj.p)) + geom_boxplot()+xlab("") +ylab("") + ylim(c(0,0.04))
gplot3 <- ggplot(data = gapdh.ljdw, aes(x="", y=peacor)) + geom_boxplot() +xlab("") +ylab("") +ylim(c(0.3,1))

gplot1.1 <- ggplot(data = try.good.ljdw.mat1, aes(x="", y=X1, fill=X1)) + geom_boxplot() +xlab("") +ylab("")+ theme(legend.position = "none")
gplot2.1 <- ggplot(data = try.good.ljdw.mat1, aes(x="", y=X2, fill=X1)) + geom_boxplot()+xlab("") +ylab("") + ylim(c(0,0.02))+ theme(legend.position = "none")
gplot3.1 <- ggplot(data = try.good.ljdw.mat1, aes(x="", y=X3, fill=X1)) + geom_boxplot() +xlab("") +ylab("") +ylim(c(0.3,1))+ theme(legend.position = "none")

gplot4 <- ggplot(data = l5dat2.na, aes(x="", y=dw.comp, fill =pcor)) + geom_boxplot() +xlab("Durbin-Watson") +ylab("") + theme(legend.position = "none")
gplot5<- ggplot(data = l5dat2.na, aes(x="", y=boxlj.p, fill = pcor)) + geom_boxplot(outlier.size=0.5)+xlab("Ljung-Box") +ylab("") + theme(legend.position = "none")+ ylim(c(0,0.25))
gplot6 <- ggplot(data = l5dat2.na, aes(x="", y=pcor, fill = pcor)) + geom_boxplot() +xlab("Pearson Cor.") +ylab("") + theme(legend.position = "none")

grid.arrange(gplot1, gplot2, gplot3, gplot1.1, gplot2.1, gplot3.1, gplot4, gplot5, gplot6, ncol=3)

#l5dat2.na[,"ind"] <- ifelse(l5dat2.na$rss <= 3000, "small", ifelse(l5dat2.na$rss >8000, "big", "ok"))
#ggplot(data = l5dat2.na, aes(x=ind, y=boxlj.p)) + geom_boxplot(aes(fill=ind), outlier.shape=NA)

load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-23a_000399.Rda")
try.good <- unlist.genparams(tst)
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-520e_001119.Rda")
try.good1 <- unlist.genparams(tst)

try.good.ljdw <- plot_sig(l5, try.good)
try.good.ljdw.mat <- data.frame(matrix(c(try.good.ljdw$dw.comp, try.good.ljdw$boxlj.p, try.good.ljdw$pearcor), ncol=3))
try.good.ljdw1 <- plot_sig(l5, try.good1)
try.good.ljdw.mat1 <- data.frame(matrix(c(try.good.ljdw1$dw.comp, try.good.ljdw1$boxlj.p, try.good.ljdw1$pearcor), ncol=3))


#miRcomp 3d
library(plotly)
plot_ly(l5dat2.na, x = ~dw.comp, y = ~pcor, z = ~boxlj.p,
             marker = list(size=1, color = ~boxlj.p, colorscale = c( '#683531','#FFE1A1'), showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Durbin-Watson'),
                      yaxis = list(title = 'Ljung-Box'),
                      zaxis = list(title = 'Pearson Correlation')))
