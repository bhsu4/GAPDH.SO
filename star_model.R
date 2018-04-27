lstar <- function(g, yc){
  1/(1+exp(-g*(yc)))
}

estar <- function(g, yc){
  (1-exp(-g*(yc)^2))
}


#plot(x=seq(-5,5,0.01), y=lstar(0, seq(-5,5,0.01)), ylim=c(0,1), type = "l", col = 1)
#lines(x=seq(-5,5,0.01), y=lstar(0.5, seq(-5,5,0.01)), ylim=c(0,1), type = "l", col = 2)
#lines(x=seq(-5,5,0.01), y=lstar(1, seq(-5,5,0.01)), type = "l", col = 3)
#lines(x=seq(-5,5,0.01), y=lstar(2.5, seq(-5,5,0.01)), type = "l", col = 4)
#lines(x=seq(-5,5,0.01), y=lstar(10, seq(-5,5,0.01)), type = "l", col = 5)
#legend("bottomright", c(paste0("Gamma", "=")), col=1:5, ncol=1, lty=1, cex=0.5)

lstar_plot <- function(xs, maxk, intk){
testk <- seq(0,maxk,intk)
plot(x=xs, y=lstar(0, xs), ylim=c(0,1), xlim=c(-2.5,2.5), type = "l", col = 1, xlab = "y_(t-d)-c", ylab = "F(y_(t-d))")
  for(i in 2:((maxk/intk)+1)){
lines(x=xs, y=lstar(testk[i], xs), ylim=c(0,1), type = "l", col = i)
  }
  legend("bottomright", c(paste0("Gamma", " =", testk[1:((maxk/intk)+1)])), col=1:((maxk/intk)+1), ncol=1, lty=1, cex=0.5)
}


star_plot <- function(xs, maxk, intk, type){
if(type == "lstar"){
    testk <- seq(0,maxk,intk)
    plot(x=xs, y=lstar(0, xs), ylim=c(0,1), xlim=c(-2.5,2.5), type = "l", col = 1, 
         xlab = expression(y[t-d]-c), ylab = expression(F(y[t-d])))
    for(i in 2:((maxk/intk)+1)){
      lines(x=xs, y=lstar(testk[i], xs), ylim=c(0,1), type = "l", col = i)
    }
    legend("bottomright", c(paste0("Gamma", " =", testk[1:((maxk/intk)+1)])), 
           col=1:((maxk/intk)+1), ncol=1, lty=1, cex=0.5, x.intersp=0.1, text.width=c(rep(0.25,5)), pt.cex = 1)
   # title(c("LSTAR With Varying Gamma"))
  }
if(type == "estar"){
  testk <- seq(0,maxk,intk)
  plot(x=xs, y=estar(0, xs), ylim=c(0,1), xlim=c(-2.5,2.5), type = "l", col = 1, 
       xlab = expression(y[t-d]-c), ylab = expression(F(y[t-d])))
  for(i in 2:((maxk/intk)+1)){
    lines(x=xs, y=estar(testk[i], xs), ylim=c(0,1), type = "l", col = i)
  }
  legend("bottomright", c(paste0("Gamma", " =", testk[1:((maxk/intk)+1)])), 
         col=1:((maxk/intk)+1), ncol=1, lty=1, cex=0.5, x.intersp=0.1, text.width=c(rep(0.25,5)), pt.cex = 1)
 # title(c("ESTAR With Varying Gamma"))
}
}

star_plot(seq(-5,5,0.01), 10, 2.5, "lstar")
star_plot(seq(-5,5,0.01), 10, 2.5, "estar")



  
  
  

