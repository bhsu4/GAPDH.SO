resid_rep_8top <- function(params, mod, rep){
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1 & as.character(mod)[[10]] == "l5"){ plot(y=resids[[rep]][1:40], 
                                                       x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                       ylim=c(-30000,30000), type="l", 
                                                       xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
      axis(2,at=seq(-25000,25000,12500)) # add a new x-axis 
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k == 1 & as.character(mod)[[10]] == "b5"){ plot(y=resids[[rep]][1:40], 
                                                            x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                            ylim=c(-30000,30000), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "l4"){ plot(y=resids[[rep]][1:40], 
                                                            x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                            ylim=c(-30000,30000), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", xaxt = "n", yaxt="n")
      #axis(2,at=seq(-25000,25000,12500)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "b4"){ plot(y=resids[[rep]][1:40], 
                                                            x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                            ylim=c(-30000,30000), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", xaxt = "n", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k > 1){
      ind2 <- 40*k
      ind1 <- ind2-39
      lines(y=resids[[rep]][ind1:ind2], x=params$fits[[rep]]$DATA$Cycles[ind1:ind2], col=k)
      abline(v=df.res1[1,], lty='dashed')
      
      #title(main= paste(names(params$fits[rep]), sub=params$fits$A$MODEL$name, sep = ", "))
      
    }
    
    else {print("no run recorded")}
  }
}

resid_rep_8bot <- function(params, mod, rep){
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1 & as.character(mod)[[10]] == "l5"){ plot(y=resids[[rep]][1:40], 
                                                       x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                       ylim=c(-30000,30000), type="l", 
                                                       xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      axis(2,at=seq(-25000,25000,12500)) # add a new x-axis 
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k == 1 & as.character(mod)[[10]] == "b5"){ plot(y=resids[[rep]][1:40], 
                                                            x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                            ylim=c(-30000,30000), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "l4"){ plot(y=resids[[rep]][1:40], 
                                                            x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                            ylim=c(-30000,30000), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      #axis(2,at=seq(-25000,25000,12500)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "b4"){ plot(y=resids[[rep]][1:40], 
                                                            x=params$fits[[1]]$DATA$Cycles[1:40], 
                                                            ylim=c(-30000,30000), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k > 1){
      ind2 <- 40*k
      ind1 <- ind2-39
      lines(y=resids[[rep]][ind1:ind2], x=params$fits[[rep]]$DATA$Cycles[ind1:ind2], col=k)
      abline(v=df.res1[1,], lty='dashed')
    }
    
    else {print("no run recorded")}
  }
}

subplot_resid_4x <- function(params){
  
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1) plot(y=resids[[k]][1:40], 
                    x=params$fits[[1]]$DATA$Cycles[1:40], 
                    type="l", ylim = c(-11000,9000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
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
                    type="l", ylim = c(-11000,9000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
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
                    type="l", ylim = c(-11000,9000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
                    xlab="", ylab="", yaxt = "n", xaxt="n")
    axis(side=2, at=seq(-8000, 8000, by=4000)) #-8000,8000 and -2000,2000
    
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
                    type="l", ylim = c(-11000,9000), # ylim=range(resids[[k]]), -3000,3000, -11000,9000
                    xlab="", ylab="", yaxt = "n", xaxt="n")
    #axis(side=2, at=seq(-8000, 8000, by=4000)) #-8000,8000 and -2000,2000
    
    if(k > 1){
      ind1 = 1 ; ind2 = 40
      #ind2 <- 40*k
      #ind1 <- ind2-39*(k)-(k-1)
      lines(y=resids[[k]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
}

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

resid_rep_8top_miR <- function(params, mod, rep){
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1 & as.character(mod)[[10]] == "l5"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                       x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                       ylim=c(-300,300), type="l", 
                                                       xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
      axis(2,at=seq(-300,300,150)) # add a new x-axis 
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k == 1 & as.character(mod)[[10]] == "b5"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                            x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                            ylim=c(-300,300), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", xaxt="n", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "l4"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                            x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                            ylim=c(-300,300), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", xaxt = "n", yaxt="n")
      #axis(2,at=seq(-25000,25000,12500)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "b4"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                            x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                            ylim=c(-300,300), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", xaxt = "n", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k > 1){
      ind2 <- max(params$fits[[rep]]$DATA$Cycles)*k
      ind1 <- ind2-(max(params$fits[[rep]]$DATA$Cycles)-1)
      lines(y=resids[[rep]][ind1:ind2], x=params$fits[[rep]]$DATA$Cycles[ind1:ind2], col=k)
      abline(v=df.res1[1,], lty='dashed')
      
      #title(main= paste(names(params$fits[rep]), sub=params$fits$A$MODEL$name, sep = ", "))
      
    }
    
    else {print("no run recorded")}
  }
}

resid_rep_8bot_miR <- function(params, mod, rep){
  for(k in 1:12){
    resids <- lapply(params$fits, resid)
    
    if(k == 1 & as.character(mod)[[10]] == "l5"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                       x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                       ylim=c(-300,300), type="l", 
                                                       xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      axis(2,at=seq(-300,300,150)) # add a new x-axis 
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k == 1 & as.character(mod)[[10]] == "b5"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                            x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                            ylim=c(-300,300), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "l4"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                            x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                            ylim=c(-300,300), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      #axis(2,at=seq(-25000,25000,12500)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    else if(k == 1 & as.character(mod)[[10]] == "b4"){ plot(y=resids[[rep]][1:max(params$fits[[rep]]$DATA$Cycles)], 
                                                            x=1:max(params$fits[[rep]]$DATA$Cycles), 
                                                            ylim=c(-300,300), type="l", 
                                                            xlab="Cycle", ylab="Fluoresence Residual", yaxt="n")
      #axis(2,at=seq(-30000,30000,15000)) # add a new x-axis
      #result ct
      df.ml1 <- modlist(params$fits[[rep]], model=mod)
      df.res1 <- getPar(df.ml1, type = "curve", cp = "cpD2", eff = "sliwin")
    }
    
    else if(k > 1){
      ind2 <- max(params$fits[[rep]]$DATA$Cycles)*k
      ind1 <- ind2-(max(params$fits[[rep]]$DATA$Cycles)-1)
      lines(y=resids[[rep]][ind1:ind2], x=params$fits[[rep]]$DATA$Cycles[ind1:ind2], col=k)
      abline(v=df.res1[1,], lty='dashed')
    }
    
    else {print("no run recorded")}
  }
}

subplot_resid_miR <- function(tester){
  
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
  
  #par(mfrow=c(2,2))
  #par(oma=c(4,4,0.5,4),mar=c(0.25,0.25,0,0))
  
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
