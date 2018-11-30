plot_sig <- function(est, listdf, macro=0, z=NULL, plot=FALSE){

  if(macro == 0 & (is.null(z)) == "TRUE"){
    par <- list() ; resids <- list() ; rss <- list()
    reg.amp <- list() ; reg.res <- list() 
    dw.amp <- list() ; dw.res <- replicate(10, list())  
    ml1 <- list() ; res1 <- list() ; paramest <- list()
    model.dw <- vector('list', 10) ; model.boxlj <- vector('list', 10) ; peacor <- vector('list', 10)
    for(i in 1:length(listdf)){
      #each replicate fit
      par[[i]] <- sub_genparams(est, listdf[[i]])
      #finding the residuals
      #for(j in 1:4){             #double loop
      #  resids[[i]][[j]] <- tryCatch({
      #    resid(par[[i]]$fits[[j]])
      #  }, error = function(e){
      #      return(replicate(max(listdf[[i]]$Cycle), NA))
      #  }) 
      #}
      try_resid <- function(x) tryCatch({resid(x)}, 
                                        error = function(e) rep(NA, max(listdf[[i]]$Cycle)))    
      resids[[i]] <- lapply(par[[i]]$fits, try_resid)
      #finding the RSS
      rss[[i]] <- lapply(resids[[i]], function(x) sum(x^2))
      #finding the Durbin-Watson stat
      reg.amp[[i]] <- lapply(listdf[[i]][ ,2:length(listdf[[i]])], 
                             function(y) dynlm(y ~ listdf[[i]][["Cycle"]])) #all amp can run
      reg.res[[i]] <- lapply(resids[[i]], function(y) try(dynlm(y ~ listdf[[i]][["Cycle"]]), 
                                                          silent=TRUE))
      #residuals contain NA
      dw.amp[[i]] <- lapply(reg.amp[[i]], function(x) durbinWatsonTest(x))
      
      for(j in 1:(length(listdf[[i]])-1)){ #NA resids cannot run
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
      
      #finding model DW
      for(j in 1:(length(listdf[[i]])-1)){ #NA resids cannot run
        model.dw[[i]][[j]] <- sum((resids[[i]][[j]]-Lag(resids[[i]][[j]], 1))^2, na.rm=TRUE)/sum(resids[[i]][[j]]^2)
        model.boxlj[[i]][[j]] <- tryCatch({Box.test(resids[[i]][[j]], lag=1, type='Ljung-Box')}, 
                 error = function(e) {return(list(statistic = NA, p.value = NA))})        
        peacor[[i]][[j]] <- cor(resids[[i]][[j]][2:max(listdf[[i]]$Cycle)], 
                                resids[[i]][[j]][1:(max(listdf[[i]]$Cycle)-1)])
      }
    }
    
    for(i in 1:length(listdf)){
      for(j in 1:(length(listdf[[i]])-1)){
        if(i==1 & j ==1){
          values <- data.frame(t(paramest[[i]][j,]))
          me <- c(dw.amp[[i]][[j]]$r, dw.amp[[i]][[j]]$dw, dw.amp[[i]][[j]]$p, dw.res[[i]][[j]]$r, 
                  dw.res[[i]][[j]]$dw, dw.res[[i]][[j]]$p, rss[[i]][[j]], res1[[i]][,j][1], res1[[i]][,j][2],
                  model.dw[[i]][[j]], model.boxlj[[i]][[j]]$statistic, model.boxlj[[i]][[j]]$p.value, peacor[[i]][[j]])
          #values[k, (length(values[k,])+1):(length(values[k,])+9)] <- cbind(values, me) 
          values <- cbind(values, t(me))
        }
        if(i==1 & j >1){
          values[j,] <- data.frame(t(paramest[[i]][j,])) #data.frame(apply(par.tst[[i]]$params[j,], c(1,2), as.numeric))
          me <- c(dw.amp[[i]][[j]]$r, dw.amp[[i]][[j]]$dw, dw.amp[[i]][[j]]$p, dw.res[[i]][[j]]$r, 
                  dw.res[[i]][[j]]$dw, dw.res[[i]][[j]]$p, rss[[i]][[j]], res1[[i]][,j][1], res1[[i]][,j][2],
                  model.dw[[i]][[j]], model.boxlj[[i]][[j]]$statistic, model.boxlj[[i]][[j]]$p.value, peacor[[i]][[j]])
          values[j, (ncol(paramest[[i]])+1):length(values)] <- t(me)
        }
        if(i > 1){
          ind2 <- i*(length(listdf[[i]])-1)
          ind1 <- ind2-(length(listdf[[i]])-2)
          values[(ind1:ind2)[j], ] <- data.frame(t(paramest[[i]][j,])) #data.frame(apply(par.tst[[i]]$params[j,], c(1,2), as.numeric))
          me <- c(dw.amp[[i]][[j]]$r, dw.amp[[i]][[j]]$dw, dw.amp[[i]][[j]]$p, 
                  dw.res[[i]][[j]]$r, dw.res[[i]][[j]]$dw, dw.res[[i]][[j]]$p, 
                  rss[[i]][[j]], res1[[i]][,j][1], res1[[i]][,j][2],
                  model.dw[[i]][[j]], model.boxlj[[i]][[j]]$statistic, 
                  model.boxlj[[i]][[j]]$p.value, peacor[[i]][[j]])
          values[(ind1:ind2)[j], (ncol(paramest[[i]])+1):length(values)] <- t(me)
        }
      }
    }
    if( est$name == "b5" || est$name == "l5"){
      names(values) <- c(c("b", "c", "d", "e", "f"), paste0(c("r", "dw", "p"), "-amp"),
                         paste0(c("r", "dw", "p"), "-res"), c("rss", "ct", "eff"), 
                         c("dw.comp", "boxlj", "boxlj.p", "pearcor"))
    }
    else{  #if( est == "l4" || est == "b4"){
      names(values) <- c(c("b", "c", "d", "e"), paste0(c("r", "dw", "p"), "-amp"),
                         paste0(c("r", "dw", "p"), "-res"), c("rss", "ct", "eff"),
                         c("dw.comp", "boxlj", "boxlj.p", "pearcor"))
    }
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
    try(points(x=res1[[i]][,k][1], y=l4_model(res1[[i]][,k][1], b=par[[i]]$params$b[k], 
                                          c=par[[i]]$params$c[k], d=par[[i]]$params$d[k], 
                                          e=par[[i]]$params$e[k]), cex=0.8, pch=16)) #CT point
 #boundaries for vertical lines
  if( (is.na(res1[[i]][,k][1]) == "TRUE") || (res1[[i]][,k][1] <= 2) || (res1[[i]][,k][1] > (max(listdf[[i]]$Cycle) - 2))){
    print(paste0(LETTERS[i], k, " ", "no ct"))
  }
  else{ 
 #small box clipping boundaries
 #   clip(min(xs), max(xs), 
 #         l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
 #                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
 #         l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
 #                                      d=par[[i]]$params$d[k], e=par[[i]]$params$e[k])) 
 #    abline(v=res1[[i]][,k][1]-2, lty=1, col=k) ; abline(v=res1[[i]][,k][1]+2, lty=1, col=k) 
 #boundaries for horizontal lines
 #   clip(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, min(listdf[[i]][ ,2:length(listdf[[i]])][[k]]),
 #                                                max(listdf[[i]][ ,2:length(listdf[[i]])][[k]])) 
 #    abline(h=l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
 #                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
 #          lty=1, col=k)
 #    abline(h=l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
 #                                          d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
 #          lty=1, col=k)
    #filling in the +/- squares for CT value
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]-2, res1[[i]][,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
    #specific small box
          #  y = c(l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
          #                                     d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
          #        l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
          #                                     d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
          #        l4_model(res1[[i]][,k][1]+2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
          #                                     d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
          #        l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
          #                                     d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]),
          #        l4_model(res1[[i]][,k][1]-2, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
          #                                     d=par[[i]]$params$d[k], e=par[[i]]$params$e[k])), 
    col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
  }
}
  if(est$name == "l5"){
    try(plot(x=xs, y=l5_model(xs, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                  d=par[[i]]$params$d[k], e=par[[i]]$params$e[k], 
                                  f=par[[i]]$params$f[k]), 
         type="l", xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf[[i]][,k+1]))), col = k, xaxt = "n"))
    points(x=xs, y=listdf[[i]][,k+1], cex=0.45) #actual points
  if(is.na(par[[i]]$params$b[k]) == "FALSE"){
    legend("topleft", c(names(listdf[[i]])[k+1]), col=k, lty=1, cex=0.65) #legend
  }
  else{}
    
    #adding box around CT values (+/- 2 cycles)
    try(points(x=res1[[i]][,k][1], y=l5_model(res1[[i]][,k][1], b=par[[i]]$params$b[k], 
                                        c=par[[i]]$params$c[k], d=par[[i]]$params$d[k], 
                                        e=par[[i]]$params$e[k], f=par[[i]]$params$f[k]), 
          cex=0.8, pch=16)) #CT point
    #boundaries for vertical lines
    if( (is.na(res1[[i]][,k][1]) == "TRUE") || (res1[[i]][,k][1] <= 2) || (res1[[i]][,k][1] > (max(listdf[[i]]$Cycle) - 2))){
      print(paste0(LETTERS[i], k, " " , "no ct"))
    }
    else{
    #filling in the +/- squares for CT value
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]-2, res1[[i]][,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
    }
  }
  if(est$name == "b4"){
    try(plot(x=xs, y=b4_model(xs, b=par[[i]]$params$b[k], c=par[[i]]$params$c[k],
                                  d=par[[i]]$params$d[k], e=par[[i]]$params$e[k]), 
         type="l", xlab="Cycle", ylab="Fluorescence", 
         ylim=c(range(unlist(listdf[[i]][,k+1]))), col = k, xaxt = "n"))
    points(x=xs, y=listdf[[i]][,k+1], cex=0.45) #actual points
  if(is.na(par[[i]]$params$b[k]) == "FALSE"){
    legend("topleft", c(names(listdf[[i]])[k+1]), col=k, lty=1, cex=0.65) #legend
  }
  else{}
    
    #adding box around CT values (+/- 2 cycles)
    try(points(x=res1[[i]][,k][1], y=b4_model(res1[[i]][,k][1], b=par[[i]]$params$b[k], 
                                          c=par[[i]]$params$c[k], d=par[[i]]$params$d[k], 
                                          e=par[[i]]$params$e[k]), cex=0.8, pch=16)) #CT point
    #boundaries for vertical lines
  if(is.na(res1[[i]][,k][1]) == "FALSE" || (res1[[i]][,k][1] <= 2) || (res1[[i]][,k][1] > (max(listdf[[i]]$Cycle) - 2))){
    #filling in the +/- squares for CT value
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]-2, res1[[i]][,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
    }
    else{
      print(paste0(LETTERS[i], k, " " , "no ct"))
    }  
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
  if(is.na(res1[[i]][,k][1]) == "FALSE" || (res1[[i]][,k][1] <= 2) || (res1[[i]][,k][1] > (max(listdf[[i]]$Cycle) - 2))){
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]-2, res1[[i]][,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
  }
    else{
      print(paste0(LETTERS[i], k, " ", "no ct"))
    }
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
      (res1[[i]][,k][1] > 38 & max(listdf[[i]]$Cycle) == 40) ||
      (res1[[i]][,k][1] > 44 & max(listdf[[i]]$Cycle) == 46)){
        print("check")
  }
  else{
 #small box
 #   clip(min(xs), max(xs), 
 #         min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), 
 #         max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]))
 #    abline(v=res1[[i]][,k][1]-2, lty=1, col=1) ; abline(v=res1[[i]][,k][1]+2, lty=1, col=1)
  #boundaries for vertical lines
 #   clip(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, min(resids[[i]][[k]]), max(resids[[i]][[k]]))
 #    abline(h=min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), lty=1, col=1)
 #    abline(h=max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), lty=1, col=1)
  #filling in the +/- squares for CT value for resids
    polygon(x = c(res1[[i]][,k][1]-2, res1[[i]][,k][1]+2, res1[[i]][,k][1]+2, 
                  res1[[i]][,k][1]-2, res1[[i]][,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
  #specific small box
 #            y = c(min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]), 
 #                  max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]),
 #                 max(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]),
 #                  min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)]),
 #                 min(resids[[i]][[k]][round(res1[[i]][,k][1]-2, digits=0):round(res1[[i]][,k][1]+2, digits=0)])),
            col= rgb(0,0,0,alpha=0.15))
        }
      } # if plot
    } #secondary list
  } #primary list
  return(values)
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
    if(is.na(par$params$b[[w]]) == "FALSE"){ #legend if non-NA's
      legend("topleft", c(names(listdf)[w+1]), col=1, lty=1, cex=0.65)
    }
      else{} #no legend if NA b/c only points show
  #adding box around CT values (+/- 2 cycles)
      try(points(x=res1[,w][1], y=l4_model(res1[,w][1], b=par$params$b[w], 
                                       c=par$params$c[w], d=par$params$d[w], 
                                       e=par$params$e[w]), cex=0.8, pch=16)) #CT point
  #boundaries for vertical lines
    if( (is.na(res1[,w][[1]]) == "TRUE") || (res1[,w][[1]] <=2) || (res1[,w][[1]] > (max(listdf$Cycle) - 2))){
      print(paste0(LETTERS[z], w, " " , "no ct")) 
    } #can't draw box for those with <=2 CT or NA
    else{
  #specific small box
#      clip(min(xs), max(xs), 
#           l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
#                                   d=par$params$d[w], e=par$params$e[w]),
#           l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
#                                   d=par$params$d[w], e=par$params$e[w])) 
#      abline(v=res1[,w][1]-2, lty=1, col=w) ; abline(v=res1[,w][1]+2, lty=1, col=w) 
  #boundaries for horizontal lines
#      clip(res1[,w][1]-2, res1[,w][1]+2, min(listdf[ ,2:length(listdf)][[w]]),
#                                         max(listdf[ ,2:length(listdf)][[w]])) 
#      abline(h=l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
#                                       d=par$params$d[w], e=par$params$e[w]), lty=1, col=w)
#      abline(h=l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
#                                       d=par$params$d[w], e=par$params$e[w]), lty=1, col=w)
  #filling in the +/- squares for CT value
      polygon(x = c(res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2, res1[,w][1]-2), 
              y = c(min(par("usr")), min(par("usr")), 
                    max(par("usr")), max(par("usr")), min(par("usr"))),
           #   y = c(l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
           #                   d=par$params$d[w], e=par$params$e[w]), 
           #          l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
           #                   d=par$params$d[w], e=par$params$e[w]),
           #          l4_model(res1[,w][1]+2, b=par$params$b[w], c=par$params$c[w],
           #                   d=par$params$d[w], e=par$params$e[w]), 
           #          l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
           #                   d=par$params$d[w], e=par$params$e[w]),
           #          l4_model(res1[,w][1]-2, b=par$params$b[w], c=par$params$c[w],
           #                   d=par$params$d[w], e=par$params$e[w])), 
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
    if(is.na(par$params$b[[w]]) == "FALSE"){ #legend if non-NA's
      legend("topleft", c(names(listdf)[w+1]), col=1, lty=1, cex=0.65)
      }
    else{}
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
  if(is.na(par$params$b[[w]]) == "FALSE"){ #legend if non-NA's
      legend("topleft", c(names(listdf)[w+1]), col=1, lty=1, cex=0.65)
  }
  else{} #no legend for NAs since no curve
      
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
  if(is.na(par$params$b[[w]]) == "FALSE"){ #legend if non-NA's
      legend("topleft", c(names(listdf)[w+1]), col=1, lty=1, cex=0.65) #col=w
  }
  else{}
      
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
         xlab = "Cycle", ylab="Fluorescence Residual")
    abline(h=0) ; points(x=res1[,w][1], y=0, cex=0.8, pch=17)
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
 #   clip(min(xs), max(xs), 
 #         min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), 
 #        max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]))
 #   abline(v=res1[,w][1]-2, lty=1, col=1) ; abline(v=res1[,w][1]+2, lty=1, col=1)
 #   clip(res1[,w][1]-2, res1[,w][1]+2, min(resids[[w]]), max(resids[[w]]))
 #   abline(h=min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), lty=1, col=1)
 #   abline(h=max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), lty=1, col=1)
    #filling in the +/- squares for CT value for resids
    polygon(x = c(res1[,w][1]-2, res1[,w][1]+2, res1[,w][1]+2, res1[,w][1]-2, res1[,w][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
 #           y = c(min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]), 
 #                 max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]),
 #                 max(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]),
 #                 min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)]),
 #                 min(resids[[w]][round(res1[,w][1]-2, digits=0):round(res1[,w][1]+2, digits=0)])),
            col= rgb(0,0,0,alpha=0.15))
      }
    } #plot bracket
  return(values)
  } #macro bracket
  
#start of macro: specific z subset
#plotting amplification curve
if(macro == 0 & z > 0){
  #define w for easier written numbers
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
  rss <- lapply(resids, function(x) sum(x^2))
  
  #finding the DW-stat
  reg.amp <- dynlm(listdf[,w+1] ~ listdf$Cycle)
  reg.res <- try(dynlm(resids[[w]] ~ listdf$Cycle), silent = TRUE)
  #finding the Durbin-Watson stat
  reg.amp <- lapply(listdf[ ,2:length(listdf)], 
                    function(y) dynlm(y ~ listdf[["Cycle"]])) #all amp can run
  reg.res <- lapply(resids, function(y) try(dynlm(y ~ listdf["Cycle"]), silent=TRUE))
  #residuals contain NA
  dw.amp <- lapply(reg.amp, function(x) durbinWatsonTest(x))
  dw.res <- list()

for(j in 1:(length(listdf)-1)){ #NA resids cannot run
  dw.res[[j]] <- tryCatch({
  durbinWatsonTest(reg.res[[j]])
  }, error=function(e) {
    return(list(r=NA, dw=NA, p=NA))
    })
  } #replace errors with NA for r, dw, p
    
  #finding the CT value
  ml1 <- modlist(listdf, model = est)
  res1 <- getPar(ml1, type = "curve", cp = "cpD2", eff = "sliwin")
  #parameter estimates
  paramest <- apply(par$params, c(1,2), as.numeric) #apply(x[k,], c(1,2), as.numeric)
    
  #all together
for(j in 1:(length(listdf)-1)){
  if(j==1){
    values <- data.frame(t(paramest[j,]))
    me <- c(dw.amp[[j]]$r, dw.amp[[j]]$dw, dw.amp[[j]]$p, dw.res[[j]]$r, 
            dw.res[[j]]$dw, dw.res[[j]]$p, rss[[j]], res1[,j][1], res1[,j][2])
  #values[k, (length(values[k,])+1):(length(values[k,])+9)] <- cbind(values, me) 
    values <- cbind(values, t(me))
    }
  if(j >1){
    values[j,] <- data.frame(t(paramest[j,])) #data.frame(apply(par.tst[[i]]$params[j,], c(1,2), as.numeric))
    me <- c(dw.amp[[j]]$r, dw.amp[[j]]$dw, dw.amp[[j]]$p, dw.res[[j]]$r, 
            dw.res[[j]]$dw, dw.res[[j]]$p, rss[[j]], res1[,j][1], res1[,j][2])
    values[j, (ncol(paramest)+1):length(values)] <- t(me) #replaces repeated paramest with dw, rss
    }
  }
  names(values) <- c(names(par$params), c("r-amp", "dw-amp", "p-amp"),
                       c("r-res", "dw-res", "p-res"), c("rss", "ct", "eff"))
  rownames(values) <- c() #getting rid of arbitrary row names  

#start plotting amp curve + resid
  for(k in 1:(length(listdf)-1)){
    xs = listdf$Cycle
  if(plot){
    if(est$name == "l4"){ #model might not be able to run
      try(plot(x=xs, y=l4_model(xs, b=par$params$b[k], c=par$params$c[k],
                                    d=par$params$d[k], e=par$params$e[k]), 
               type="l", xlab="Cycle", ylab="Fluorescence", 
               ylim=c(range(unlist(listdf[,k+1]))), col = k, xaxt = "n"))
          points(x=xs, y=listdf[,k+1], cex=0.45) #actual points
          
    if(is.na(par$params$b[k]) == "FALSE"){
      legend("topleft", c(names(listdf)[k+1]), col=k, lty=1, cex=0.65) #legend
      } #adds legend for line of model
    else{} #only add legend if able to run model
          
  #adding box around CT values (+/- 2 cycles)
    try(points(x=res1[,k][1], y=l4_model(res1[,k][1], b=par$params$b[k], 
                                         c=par$params$c[k], d=par$params$d[k], 
                                         e=par$params$e[k]), cex=0.8, pch=16)) #CT point
  #no CT value for near end CTs or NAs
    if( (is.na(res1[,k][1]) == "TRUE") || (res1[,k][1] <= 2)){
         print(paste0(LETTERS[i], k, " ", "no ct"))
      }
    else{ 
    #big box around +/- 2 cycles
    polygon(x = c(res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, 
                  res1[,k][1]-2, res1[,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
          }
        }
        
  if(est$name == "l5"){ #model might not be able to run
    try(plot(x=xs, y=l5_model(xs, b=par$params$b[k], c=par$params$c[k],
                                  d=par$params$d[k], e=par$params$e[k], f=par$params$f[k]), 
             type="l", xlab="Cycle", ylab="Fluorescence", 
             ylim=c(range(unlist(listdf[,k+1]))), col = k, xaxt = "n")) #col=k
        points(x=xs, y=listdf[,k+1], cex=0.45) #actual points
  if(is.na(par$params$b[k]) == "FALSE"){
    legend("topleft", c(names(listdf)[k+1]), col=k, lty=1, cex=0.65) #legend
    } #adds legend for line of model
  else{} #only add legend if able to run model
    #adding box around CT values (+/- 2 cycles)
    try(points(x=res1[,k][1], y=l5_model(res1[,k][1], b=par$params$b[k], 
                                         c=par$params$c[k], d=par$params$d[k], 
                                         e=par$params$e[k], f=par$params$f[k]), 
        cex=0.8, pch=16)) #CT point
    #no CT value for near end CTs or NAs
  if( (is.na(res1[,k][1]) == "TRUE")  || (res1[,k][1] <= 2) || (res1[,k][1] > (max(listdf$Cycle) - 2))){
       print(paste0(LETTERS[i], k, " ", "no ct"))
      }
  else{ 
    #big box around +/- 2 cycles
    polygon(x = c(res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, 
                  res1[,k][1]-2, res1[,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
      }
    }  
        
  if(est$name == "b4"){ #model might not be able to run
      try(plot(x=xs, y=b4_model(xs, b=par$params$b[k], c=par$params$c[k],
                                    d=par$params$d[k], e=par$params$e[k]), 
          type="l", xlab="Cycle", ylab="Fluorescence", 
          ylim=c(range(unlist(listdf[,k+1]))), col = k, xaxt = "n")) #col=k
      points(x=xs, y=listdf[,k+1], cex=0.45) #actual points
  if(is.na(par$params$b[k]) == "FALSE"){
      legend("topleft", c(names(listdf)[k+1]), col=k, lty=1, cex=0.65) #legend #col=k
    } #adds legend for line of model
  else{} #only add legend if able to run model
  #adding box around CT values (+/- 2 cycles)
  try(points(x=res1[,k][1], y=b4_model(res1[,k][1], b=par$params$b[k], 
                                       c=par$params$c[k], d=par$params$d[k], 
                                       e=par$params$e[k]), cex=0.8, pch=16)) #CT point
  #no CT value for near end CTs or NAs
  if(is.na(res1[,k][1]) == "FALSE"){
    #big box around +/- 2 cycles
    polygon(x = c(res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, 
                  res1[,k][1]-2, res1[,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
      }
  else{ 
    print(paste0(LETTERS[i], k, " ", "no ct"))
      }
    }
        
  if(est$name == "b5"){ #model might not be able to run
    try(plot(x=xs, y=b5_model(xs, b=par$params$b[k], c=par$params$c[k],
                                  d=par$params$d[k], e=par$params$e[k], f=par$params$f[k]), 
        type="l", xlab="Cycle", ylab="Fluorescence", 
        ylim=c(range(unlist(listdf[,k+1]))), col = k, xaxt = "n")) #col=k
        points(x=xs, y=listdf[,k+1], cex=0.45) #actual points
          
  if(is.na(par$params$b[k]) == "FALSE"){
    legend("topleft", c(names(listdf)[k+1]), col=k, lty=1, cex=0.65) #legend #col=k
      } #adds legend for line of model
  else{} #only add legend if able to run model
          
    #adding box around CT values (+/- 2 cycles)
    try(points(x=res1[,k][1], y=b5_model(res1[,k][1], b=par$params$b[k], 
                                         c=par$params$c[k], d=par$params$d[k], 
                                         e=par$params$e[k], f=par$params$f[k]), 
        cex=0.8, pch=16)) #CT point
    #no CT value for near end CTs or NAs
  if( (is.na(res1[,k][1]) == "FALSE")){
    #big box around +/- 2 cycles
    polygon(x = c(res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, 
                  res1[,k][1]-2, res1[,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15)) #density=10, angle=-45, col = "grey", lty=2)
      }
  else{ 
    print(paste0(LETTERS[i], k, " ", "no ct"))
    }
  } #b5
        
  if(unique(is.na(resids[[k]])) == "FALSE"){
    plot(y = resids[[k]][1:length(listdf$Cycle)], 
         x = par$fits[[k]]$DATA$Cycles[1:length(listdf$Cycle)], 
         ylim = range(unlist(resids[[k]])), 
         xlab = "Cycle", ylab="Fluorescence Residual")
         abline(h=0) ; points(x=res1[,k][1], y=0, cex=0.8, pch=17)
    }
  else{
    plot(1, type="n" , xlab="n", ylab="", 
         xlim=c(0, sum(is.na(resids[[k]]))), ylim=c(0,1))
    }
  #boundaries for horizontal lines
  if( (res1[,k][1] <= 2) || (is.na(res1[,k][1])) == "TRUE" ||
      (res1[,k][1] > 38 & max(listdf$Cycle) == 40) ||
      (res1[,k][1] > 44 & max(listdf$Cycle) == 46)){
    print("check")
    }
  else{
    polygon(x = c(res1[,k][1]-2, res1[,k][1]+2, res1[,k][1]+2, 
                  res1[,k][1]-2, res1[,k][1]-2), 
            y = c(min(par("usr")), min(par("usr")), 
                  max(par("usr")), max(par("usr")), min(par("usr"))),
            col= rgb(0,0,0,alpha=0.15))
        }
      } # if plot
    } #for j in 1:(length-1) subsetted list
    return(values)
  }#macro==0 z>0
}

sig_est <- function(est, orgdata, getfiles){
  #get the files which overwrites tst 
  files <- getfiles
  targnames <- unique(files) #all targets
  targnames <- gsub(".Rda", "", targnames) #remove .Rda to match later
  sampnames <- sort(unique(c(orgdata$SampleID)))
  #repeat all sample names for the chosen targets in files
  tmp <- rep(sampnames, length(files))
  #creating result data.frame
  targlength <- length(targnames) #length of all targets
  replength <- length(unique(orgdata$SampleID)) #sum(lengths(tst)) #length of total reps for each target
  if(est$name == "l4" || est$name == "b4"){
    res <- data.frame(
      #target categories
      TargetName = rep(targnames, each = replength), SampleID = rep(sampnames, targlength), 
      Group = gsub("_." , "", tmp), FeatureSet = rep(NA, each = replength), 
      #parameter est for 4 parm
      b = rep(NA, targlength * replength), c = rep(NA, targlength * replength), 
      d = rep(NA, targlength * replength), e = rep(NA, targlength * replength),
      #dw statistics
      r.amp = rep(NA, targlength * replength), dw.amp = rep(NA, targlength * replength), p.amp = rep(NA, targlength * replength), 
      r.res = rep(NA, targlength * replength), dw.res = rep(NA, targlength * replength), p.res = rep(NA, targlength * replength), 
      #rss and getPar statistics
      rss = rep(NA, targlength * replength), ct = rep(NA, targlength * replength), eff = rep(NA,targlength * replength),
      dw.comp = rep(NA, targlength * replength), boxlj.x2 = rep(NA, targlength * replength), 
      boxlj.p = rep(NA,targlength * replength), pcor = rep(NA, targlength * replength)
    )
  }
  if(est$name == "l5" || est$name == "b5"){
    res <- data.frame(
      #target categories
      TargetName = rep(targnames, each = replength), SampleID = rep(sampnames, targlength), 
      Group = gsub("_." , "", tmp), FeatureSet = rep(NA, each = replength), 
      #parameter est for 5 parm
      b = rep(NA, targlength * replength), c = rep(NA, targlength * replength), d = rep(NA, targlength * replength), 
      e = rep(NA, targlength * replength), f = rep(NA, targlength * replength), 
      #dw statistics
      r.amp = rep(NA, targlength * replength), dw.amp = rep(NA, targlength * replength), p.amp = rep(NA, targlength * replength), 
      r.res = rep(NA, targlength * replength), dw.res = rep(NA, targlength * replength), p.res = rep(NA, targlength * replength), 
      #rss and getPar statistics
      rss = rep(NA, targlength * replength), ct = rep(NA, targlength * replength), eff = rep(NA,targlength * replength),
      dw.comp = rep(NA, targlength * replength), boxlj.x2 = rep(NA, targlength * replength), 
      boxlj.p = rep(NA,targlength * replength), pcor = rep(NA, targlength * replength)
      
    )
  }
  foreach(k=1:length(files)) %do% {    #(k in 1:length(files)){  #foreach (k = 1:length(files)) %do% {
    load(file = files[[k]])
    print(getfiles[[k]])
    try <- unlist.genparams(tst)
    ind2 <- length(unique(orgdata$SampleID))*k  ; ind1 <- ind2-(length(unique(orgdata$SampleID))-1)
    if((grepl(tst[[1]][[1]]$TargetName[1], res[,"TargetName"][ind1]) == "TRUE") & 
       (gsub("[[:alpha:]]","", est$name) == "5") == "TRUE"){
      res[ind1:ind2, 5:22] <- plot_sig(est, try)
      res[ind1:ind2, "FeatureSet"] <- rep(as.character(unique(lapply(tst[[1]], 
                                                    function(x) unique(x$FeatureSet)))), replength)
    }
    if((grepl(tst[[1]][[1]]$TargetName[1], res[,"TargetName"][ind1]) == "TRUE") & 
       (gsub("[[:alpha:]]","", est$name) == "4") == "TRUE"){
      res[ind1:ind2, 5:21] <- plot_sig(est, try)
      res[ind1:ind2, "FeatureSet"] <- rep(as.character(unique(lapply(tst[[1]], 
                                                     function(x) unique(x$FeatureSet)))), replength) 
    }
  }
  return(res)
}


