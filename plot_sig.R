plot_sig <- function(est, listdf, macro=0, z=NULL, plot=FALSE){

if(macro == 0){
  par <- list() ; resids <- replicate(10, list()); rss <- list()
  reg.amp <- list() ; reg.res <- list() 
  dw.amp <- list() ; dw.res <- replicate(10, list())  
  ml1 <- list() ; res1 <- list() ; paramest <- list()
for(i in 1:length(listdf)){
  #each replicate fit
  par[[i]] <- sub_genparams(est, listdf[[i]])
  #finding the residuals
#for(j in 1:4){ double loop
#  resids[[i]][[j]] <- tryCatch({
#    resid(par[[i]]$fits[[j]])
#  }, error = function(e){
#      return(replicate(max(listdf[[i]]$Cycle), NA))
#  }) 
#}
  try_resid <- function(x) tryCatch({resid(x)}, 
                           error = function(e) rep(NA, max(try[[i]]$Cycle)))    
  resids[[i]] <- lapply(par.tst[[i]]$fits, try_resid)
  
  #finding the RSS
  rss[[i]] <- lapply(resids[[i]], function(x) sum(x^2))
  #finding the Durbin-Watson stat
  reg.amp[[i]] <- lapply(listdf[[i]][ ,2:length(listdf[[i]])], 
                         function(y) dynlm(y ~ listdf[[i]][["Cycle"]])) #all amp can run
  reg.res[[i]] <- lapply(resids[[i]], function(y) try(dynlm(y ~ listdf[[i]][["Cycle"]]), 
                                                      silent=TRUE))
  #residuals contain NA
  dw.amp[[i]] <- lapply(reg.amp[[i]], function(x) durbinWatsonTest(x))

for(j in 1:4){ #NA resids cannot run
  dw.res[[i]][[j]] <- tryCatch({
    durbinWatsonTest(reg.res[[i]][[j]])
  }, error=function(e) {
    return(list(r=NA, dw=NA, p=NA))
  })
} #replace errors with NA for r, dw, p
  
  #finding CT value
  ml1[[i]] <- modlist(listdf[[i]], model = est)
  res1[[i]] <- getPar(ml1[[i]], type = "curve", cp = "cpD2", eff = "sliwin")
  #parameter estimates
  paramest[[i]] <- apply(par[[i]]$params, c(1,2), as.numeric) #apply(x[k,], c(1,2), as.numeric)
}
for(i in 1:length(listdf)){
  for(j in 1:(length(listdf[[i]])-1)){
    if(i==1 & j ==1){
      values <- data.frame(t(paramest[[i]][j,]))
      me <- c(dw.amp[[i]][[j]]$r, dw.amp[[i]][[j]]$dw, dw.amp[[i]][[j]]$p, dw.res[[i]][[j]]$r, 
              dw.res[[i]][[j]]$dw, dw.res[[i]][[j]]$p, rss[[i]][[j]], res1[[i]][,j][1], res1[[i]][,j][2])
      #values[k, (length(values[k,])+1):(length(values[k,])+9)] <- cbind(values, me) 
      values <- cbind(values, t(me))
    }
    if(i==1 & j >1){
      values[j,] <- data.frame(t(paramest[[i]][j,])) #data.frame(apply(par.tst[[i]]$params[j,], c(1,2), as.numeric))
      me <- c(dw.amp[[i]][[j]]$r, dw.amp[[i]][[j]]$dw, dw.amp[[i]][[j]]$p, dw.res[[i]][[j]]$r, 
              dw.res[[i]][[j]]$dw, dw.res[[i]][[j]]$p, rss[[i]][[j]], res1[[i]][,j][1], res1[[i]][,j][2])
      values[j, (ncol(paramest[[i]])+1):length(values)] <- t(me)
    }
    if(i > 1){
      ind2 <- i*(length(listdf[[i]])-1)
      ind1 <- ind2-(length(listdf[[i]])-2)
      values[(ind1:ind2)[j], ] <- data.frame(t(paramest[[i]][j,])) #data.frame(apply(par.tst[[i]]$params[j,], c(1,2), as.numeric))
      me <- c(dw.amp[[i]][[j]]$r, dw.amp[[i]][[j]]$dw, dw.amp[[i]][[j]]$p, 
              dw.res[[i]][[j]]$r, dw.res[[i]][[j]]$dw, dw.res[[i]][[j]]$p, 
              rss[[i]][[j]], res1[[i]][,j][1], res1[[i]][,j][2])
      values[(ind1:ind2)[j], (ncol(paramest[[i]])+1):length(values)] <- t(me)
    }
  }
}
names(values) <- c(c("b", "c", "d", "e", "f"), paste0(c("r", "dw", "p"), "-amp"),
                   paste0(c("r", "dw", "p"), "-res"), c("rss", "ct", "eff"))
rownames(values) <- c() #getting rid of arbitrary row names

#two graphs on top each other
par(mfrow=c(2,1))
par(oma=c(4,4,4,4),mar=c(0.25,0.25,0,0))

#plotting amplification curve
for(i in 1:length(listdf)){
  for(k in 1:(length(listdf[[i]])-1)){
    xs = listdf[[i]]$Cycle
if(plot){
  if(est$name == "l4"){ #model might not be able to run
    try(plot(x=xs, y=l4_model(xs, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                  d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
         type="l", xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf[[i]][,k+1]))), col = k, xaxt = "n"))
         points(x=xs, y=listdf[[i]][,k+1], cex=0.45) #actual points
  if(is.na(par[[i]]$params$b[k]) == "FALSE"){
         legend("topleft", c(names(listdf[[i]])[k+1]), col=k, lty=1, cex=0.65) #legend
  } #adds legend for line of model
         else{} #only add legend if able to run model
 #adding box around CT values (+/- 2 cycles)
    points(x=res1[[i]][,k][1], y=l4_model(res1[[i]][,k][1], b=par[[i]]$params$b[k], 
                                          c=par[[i]]$params$c[k], d=par[[i]]$params$d[k], 
                                          e=par[[i]]$params$e[k]), cex=0.8, pch=16) #CT point
 #boundaries for vertical lines
    clip(min(xs), max(xs), 
         l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
         l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k])) 
    abline(v=res1[[i]][,k][1]-2, lty=1, col=k) ; abline(v=res1[[i]][,k][1]+2, lty=1, col=k) 
 #boundaries for horizontal lines
    clip(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, min(listdf[[i]][ ,2:length(listdf[[i]])][[k]]),
                                                 max(listdf[[i]][ ,2:length(listdf[[i]])][[k]])) 
    abline(h=l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
           lty=1, col=k)
    abline(h=l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
           lty=1, col=k)
 #filling in the +/- squares for CT value
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]+2, res1[[i]][,k][1]-2), 
            y = c(l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
                  l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
                  l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
                  l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
                  l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k])), 
    col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
  }
  if(est$name == "l5"){
    plot(x=xs, y=l5_model(xs, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                              d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                              f=par[[i]]$params$f[k]), 
         type="l", xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf[[i]][,k+1]))), col = k, xaxt = "n")
    points(x=xs, y=listdf[[i]][,k+1], cex=0.45) #actual points
    legend("topleft", c(names(listdf[[i]])[k+1]), col=k, lty=1, cex=0.65) #legend
    
    #adding box around CT values (+/- 2 cycles)
    points(x=res1[[i]][,k][1], y=l5_model(res1[[i]][,k][1], b=par[[i]]$params$b[k], 
                                          c=par[[i]]$params$c[k], d=par[[i]]$params$d[k], 
                                          e=par[[i]]$params$e[k], f=par[[i]]$params$f[k]), 
          cex=0.8, pch=16) #CT point
    #boundaries for vertical lines
    clip(min(xs), max(xs), 
         l5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                      f=par[[i]]$params$f[k]),
         l5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                      f=par[[i]]$params$f[k])) 
    abline(v=res1[[i]][,k][1]-2, lty=1, col=k) ; abline(v=res1[[i]][,k][1]+2, lty=1, col=k) 
    #boundaries for horizontal lines
    clip(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, min(listdf[[i]][ ,2:length(listdf[[i]])][[k]]),
         max(listdf[[i]][ ,2:length(listdf[[i]])][[k]])) 
    abline(h=l5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                          f=par[[i]]$params$f[k]), lty=1, col=k)
    abline(h=l5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                          f=par[[i]]$params$f[k]), lty=1, col=k)
    #filling in the +/- squares for CT value
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]+2, res1[[i]][,k][1]-2), 
            y = c(l5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                               f=par[[i]]$params$f[k]), 
                  l5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                               f=par[[i]]$params$f[k]),
                  l5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                               f=par[[i]]$params$f[k]), 
                  l5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                               f=par[[i]]$params$f[k]),
                  l5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                               f=par[[i]]$params$f[k])), 
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
  } 
  if(est$name == "b4"){
    plot(x=xs, y=b4_model(xs, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                              d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
         type="l", xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf[[i]][,k+1]))), col = k, xaxt = "n")
    points(x=xs, y=listdf[[i]][,k+1], cex=0.45) #actual points
    legend("topleft", c(names(listdf[[i]])[k+1]), col=k, lty=1, cex=0.65) #legend
    
    #adding box around CT values (+/- 2 cycles)
    points(x=res1[[i]][,k][1], y=b4_model(res1[[i]][,k][1], b=par[[i]]$params$b[k], 
                                          c=par[[i]]$params$c[k], d=par[[i]]$params$d[k], 
                                          e=par[[i]]$params$e[k]), cex=0.8, pch=16) #CT point
    #boundaries for vertical lines
    clip(min(xs), max(xs), 
         b4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
         b4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k])) 
    abline(v=res1[[i]][,k][1]-2, lty=1, col=k) ; abline(v=res1[[i]][,k][1]+2, lty=1, col=k) 
    #boundaries for horizontal lines
    clip(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, min(listdf[[i]][ ,2:length(listdf[[i]])][[k]]),
                                                 max(listdf[[i]][ ,2:length(listdf[[i]])][[k]])) 
    abline(h=b4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
           lty=1, col=k)
    abline(h=b4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
           lty=1, col=k)
    #filling in the +/- squares for CT value
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]+2, res1[[i]][,k][1]-2), 
            y = c(b4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
                  b4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
                  b4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
                  b4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
                  b4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k])), 
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
  }
  if(est$name == "b5"){
    try(plot(x=xs, y=b5_model(xs, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                  d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                  f=par[[i]]$params$f[k]), 
         type="l", xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf[[i]][,k+1]))), col = k, xaxt = "n"))
    points(x=xs, y=listdf[[i]][,k+1], cex=0.45) #actual points
  if(is.na(par[[i]]$params$b[k]) == "FALSE"){
    legend("topleft", c(names(listdf[[i]])[k+1]), col=k, lty=1, cex=0.65) #legend
  }
  else{} #no legend if model fails to run
    
    #adding box around CT values (+/- 2 cycles)
    try(points(x=res1[[i]][,k][1], y=b5_model(res1[[i]][,k][1], b=par[[i]]$params$b[k], 
                                        c=par[[i]]$params$c[k], d=par[[i]]$params$d[k], 
                                        e=par[[i]]$params$e[k], f=par[[i]]$params$f[k]), 
           cex=0.8, pch=16)) #CT point
    #boundaries for vertical lines
  if(is.na(res1[[i]][,k][1]) == "FALSE"){
    clip(min(xs), max(xs), 
         b5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                      f=par[[i]]$params$f[k]),
         b5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                      f=par[[i]]$params$f[k])) 
    abline(v=res1[[i]][,k][1]-2, lty=1, col=k) ; abline(v=res1[[i]][,k][1]+2, lty=1, col=k) 
    #boundaries for horizontal lines
    clip(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, min(listdf[[i]][ ,2:length(listdf[[i]])][[k]]),
         max(listdf[[i]][ ,2:length(listdf[[i]])][[k]])) 
    abline(h=b5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                          f=par[[i]]$params$f[k]), lty=1, col=k)
    abline(h=b5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                          f=par[[i]]$params$f[k]), lty=1, col=k)
    #filling in the +/- squares for CT value
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]+2, res1[[i]][,k][1]-2), 
            y = c(b5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                               f=par[[i]]$params$f[k]), 
                  b5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                               f=par[[i]]$params$f[k]),
                  b5_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                               f=par[[i]]$params$f[k]), 
                  b5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                               f=par[[i]]$params$f[k]),
                  b5_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                               d=par[[i]]$params$d[k], e=par[[i]]$params$e[k],
                                               f=par[[i]]$params$f[k])), 
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
  }
    else{
      print(paste0(LETTERS[i], k, " ", "no ct"))
    }
  #plotting residuals
  if(unique(is.na(resids[[i]][[k]])) == "FALSE"){
    plot(y = resids[[i]][[k]][1:length(listdf[[i]]$Cycle)], 
         x = par[[i]]$fits[[k]]$DATA$Cycles[1:length(listdf[[i]]$Cycle)], 
         ylim = range(unlist(resids[[i]][[k]])), 
         xlab = "Cycle", ylab="Fluorescence Residual")
         abline(h=0) ; points(x=res1[[i]][,k][1], y=0, cex=0.8, pch=17)
  }
  else{
    plot(1, type="n" , xlab="n", ylab="", 
         xlim=c(0, sum(is.na(resids[[i]][[k]]))), ylim=c(0,1))
  }
  #boundaries for horizontal lines
  if( (res1[[i]][,k][1] <= 2) || (is.na(res1[[i]][,k][1])) == "TRUE" ||
      (res1[[i]][,k][1] > 38 & max(try[[i]]$Cycle) == 40) ||
      (res1[[i]][,k][1] > 44 & max(try[[i]]$Cycle) == 46)){
        print("check")
  }
  else{
    clip(min(xs), max(xs), 
         min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), 
         max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]))
    abline(v=res1[[i]][,k][1]-2, lty=1, col=1) ; abline(v=res1[[i]][,k][1]+2, lty=1, col=1)
  #boundaries for vertical lines
    clip(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, min(resids[[i]][[k]]), max(resids[[i]][[k]]))
    abline(h=min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), lty=1, col=1)
    abline(h=max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), lty=1, col=1)
  #filling in the +/- squares for CT value for resids
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, res1[[i]][,k][1]+2, res1[[i]][,k][1]-2), 
            y = c(min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), 
                  max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]),
                  max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]),
                  min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]),
                  min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)])),
            col= rgb(0,0,0,alpha=0.15))
          }
        } #if b5 bracket
      } # if plot
    } #secondary list
  } #primary list
} #macro bracket
  
#start of microscoping: specific w term evaluation

#plotting amplification curve
if(macro > 0){
  w = macro
  #re-state specified z-list
  listdf <- listdf[[z]]
  #getting parameter estimates
  par <- sub_genparams(est, listdf)
  #finding the residuals
  resids <- lapply(par$fits, resid)
  #finding the RSS
  rss <- sum(resids[[w]]^2)
  #finding the DW-stat
  reg.amp <- dynlm(listdf[,w+1] ~ listdf$Cycle)
  reg.res <- dynlm(resids[[w]] ~ listdf$Cycle)
  dw.amp <- durbinWatsonTest(reg.amp) ; dw.res <- durbinWatsonTest(reg.res)
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
    if(est$name == "l4"){
      plot(x=xs, y=l4_model(xs, b=par$params$b[w], c=par$params$c[w],
                                d=par$params$d[w], e=par$params$e[w]), type="l",  
           xlab="Cycle", ylab="Fluorescence", 
           ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                        names(unlist(listdf)))])])), col = w, xaxt = "n")
      points(x=xs, y=listdf[,w+1], cex=0.45)
      legend("topleft", c(names(listdf)[w+1]), col=w, lty=1, cex=0.65)
  #adding box around CT values (+/- 2 cycles)
      points(x=res1[,w][1], y=l4_model(res1[,w][1], b=par$params$b[w], 
                                       c=par$params$c[w], d=par$params$d[w], 
                                       e=par$params$e[w]), cex=0.8, pch=16) #CT point
  #boundaries for vertical lines
      clip(min(xs), max(xs), 
           l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w]),
           l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w])) 
      abline(v=res1[,w][1]-2, lty=1, col=w) ; abline(v=res1[,w][1]+2, lty=1, col=w) 
  #boundaries for horizontal lines
      clip(res1[,w][1]-2, res1[,w][1]+2, min(listdf[ ,2:length(listdf)][[w]]),
                                         max(listdf[ ,2:length(listdf)][[w]])) 
      abline(h=l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w]), lty=1, col=w)
      abline(h=l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w]), lty=1, col=w)
  #filling in the +/- squares for CT value
      polygon(x = c(res1[,w][1]-2, res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2), 
              y = c(l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w]), 
                    l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w]),
                    l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w]), 
                    l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w]),
                    l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w])), 
              col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
    }
    if(est$name == "b4"){
      plot(x=xs, y=b4_model(xs, b=par$params$b[w], c=par$params$c[w],
                                d=par$params$d[w], e=par$params$e[w]), type="l",  
           xlab="Cycle", ylab="Fluorescence", 
           ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                        names(unlist(listdf)))])])), col = w, xaxt = "n")
      points(x=xs, y=listdf[,w+1], cex=0.45)
      legend("topleft", c(names(listdf)[w+1]), col=w, lty=1, cex=0.65)
  #adding box around CT values (+/- 2 cycles)
      points(x=res1[,w][1], y=b4_model(res1[,w][1], b=par$params$b[w], 
                                       c=par$params$c[w], d=par$params$d[w], 
                                       e=par$params$e[w]), cex=0.8, pch=16) #CT point
  #boundaries for vertical lines
      clip(min(xs), max(xs), 
           b4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w]),
           b4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w])) 
      abline(v=res1[,w][1]-2, lty=1, col=w) ; abline(v=res1[,w][1]+2, lty=1, col=w) 
  #boundaries for horizontal lines
      clip(res1[,w][1]-2, res1[,w][1]+2, min(listdf[ ,2:length(listdf)][[w]]),
                                         max(listdf[ ,2:length(listdf)][[w]])) 
      abline(h=b4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w]), lty=1, col=w)
      abline(h=b4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w]), lty=1, col=w)
  #filling in the +/- squares for CT value
      polygon(x = c(res1[,w][1]-2, res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2), 
              y = c(b4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w]), 
                    b4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w]),
                    b4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w]), 
                    b4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w]),
                    b4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w])), 
              col= rgb(0,0,0,alpha=0.15))
    }
  if(est$name == "l5"){
      plot(x=xs, y=l5_model(xs, b=par$params$b[w], c=par$params$c[w],
                                d=par$params$d[w], e=par$params$e[w],
                                f=par$params$f[w]), type="l",  
           xlab="Cycle", ylab="Fluorescence", 
           ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                        names(unlist(listdf)))])])), col = w, xaxt = "n")
      points(x=xs, y=listdf[,w+1], cex=0.45)
      legend("topleft", c(names(listdf)[w+1]), col=w, lty=1, cex=0.65)
  #adding box around CT values (+/- 2 cycles)
      points(x=res1[,w][1], y=l5_model(res1[,w][1], b=par$params$b[w], #CT point
                                       c=par$params$c[w], d=par$params$d[w], 
                                       e=par$params$e[w], f=par$params$f[w]), cex=0.8, pch=16) 
  #boundaries for vertical lines
      clip(min(xs), max(xs), 
           l5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w], f=par$params$f[w]),
           l5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w], f=par$params$f[w])) 
      abline(v=res1[,w][1]-2, lty=1, col=46) ; abline(v=res1[,w][1]+2, lty=1, col=w) 
  #boundaries for horizontal lines
      clip(res1[,w][1]-2, res1[,w][1]+2, min(listdf[ ,2:length(listdf)][[w]]),
                                         max(listdf[ ,2:length(listdf)][[w]])) 
      abline(h=l5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w], 
                                       f=par$params$f[w]), lty=1, col=w)
      abline(h=l5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w],
                                       f=par$params$f[w]), lty=1, col=w)
  #filling in the +/- squares for CT value
      polygon(x = c(res1[,w][1]-2, res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2), 
              y = c(l5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w],
                                            f=par$params$f[w]), 
                    l5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w],
                                            f=par$params$f[w]),
                    l5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w],
                                            f=par$params$f[w]), 
                    l5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w], 
                                            f=par$params$f[w]),
                    l5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                            d=par$params$d[w], e=par$params$e[w], 
                                            f=par$params$f[w])), 
              col= rgb(0,0,0,alpha=0.15))
    }
  if(est$name == "b5"){ 
      plot(x=xs, y=b5_model(xs, b=par$params$b[w], c=par$params$c[w],
                                d=par$params$d[w], e=par$params$e[w],
                                f=par$params$f[w]), type="l",  
           xlab="Cycle", ylab="Fluorescence", 
           ylim=c(range(unlist(listdf)[(names(unlist(listdf))[!grepl("Cycle", 
                                        names(unlist(listdf)))])])), col = w, xaxt = "n")
      points(x=xs, y=listdf[,w+1], cex=0.45)
      legend("topleft", c(names(listdf)[w+1]), col=w, lty=1, cex=0.65)
   #adding box around CT values (+/- 2 cycles)
      points(x=res1[,w][1], y=b5_model(res1[,w][1], b=par$params$b[w], #CT point
                                       c=par$params$c[w], d=par$params$d[w], 
                                       e=par$params$e[w], f=par$params$f[w]), cex=0.8, pch=16) 
   #boundaries for vertical lines
      clip(min(xs), max(xs), 
           b5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w],
                                   f=par$params$f[w]),
           b5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                   d=par$params$d[w], e=par$params$e[w],
                                   f=par$params$f[w])) 
      abline(v=res1[,w][1]-2, lty=1, col=46) ; abline(v=res1[,w][1]+2, lty=1, col=w) 
   #boundaries for horizontal lines
      clip(res1[,w][1]-2, res1[,w][1]+2, min(listdf[ ,2:length(listdf)][[w]]),
                                         max(listdf[ ,2:length(listdf)][[w]])) 
      abline(h=b5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w],
                                       f=par$params$f[w]), lty=1, col=w)
      abline(h=b5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                                       d=par$params$d[w], e=par$params$e[w],
                                       f=par$params$f[w]), lty=1, col=w) 
   #filling in the +/- squares for CT value
      polygon(x = c(res1[,w][1]-2, res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2), 
              y = c(b5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w], f=par$params$f[w]), 
                    b5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w], f=par$params$f[w]),
                    b5_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w], f=par$params$f[w]), 
                    b5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w], f=par$params$f[w]),
                    b5_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
                             d=par$params$d[w], e=par$params$e[w], f=par$params$f[w])), 
              col= rgb(0,0,0,alpha=0.15))
    }
  #plotting residuals
    plot(y = resids[[w]][1:length(listdf$Cycle)], 
         x = par$fits[[w]]$DATA$Cycles[1:length(listdf$Cycle)], 
         ylim = range(unlist(resids)), 
         xlab = "Cycle", ylab="Fluorescence Residual")
    abline(h=0) ; points(x=res1[,w][1], y=0, cex=0.8, pch=17)
    clip(min(xs), max(xs), 
         min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), 
         max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]))
    abline(v=res1[,w][1]-2, lty=1, col=1) ; abline(v=res1[,w][1]+2, lty=1, col=1)
    clip(res1[,w][1]-2, res1[,w][1]+2, min(resids[[w]]), max(resids[[w]]))
    abline(h=min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), lty=1, col=1)
    abline(h=max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), lty=1, col=1)
    #filling in the +/- squares for CT value for resids
    polygon(x = c(res1[,w][1]-2, res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2), 
            y = c(min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), 
                  max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]),
                  max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]),
                  min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]),
                  min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)])),
            col= rgb(0,0,0,alpha=0.15))
    }
  return(values)
  }
}


