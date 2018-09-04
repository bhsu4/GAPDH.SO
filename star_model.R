#keeps first column cycle, log10 rest of fluorescence 
library(tidyverse) ; library(gtools) ; library(strucchange) ; library(data.table)

hello <- function(subs){

subslog <- lapply(subs, function(x) x %>% 
                  mutate_each(funs(log10(.)), 2:13)) #original log10 terms 
subslogcb <- lapply(subslog, function(x) x %>% 
                      select(-starts_with("Cycle")) %>% #delete cycle lag
                      mutate_each(funs(lag(., k=2)), 1:12))  #lagged terms
subslogcb <- lapply(subslogcb, function(x) setNames(x, paste0(colnames(x), "_lag"))) #colnames for lagged
subs.all <- Map(cbind, subslog, subslogcb) #mapped for column bind for log and logcb
subs.all <- lapply(subs.all, function(x) x %>% 
                     select(noquote(mixedsort(colnames(.)))) %>% #correct order except for Cycles
                     select(starts_with("Cycle"), everything())) #correct order
#subs.all <- lapply(subs.all, function(x) x %>% select(-starts_with("Cycle_")))
brks <- function(subsfcn){ #function for subset breakdate est.
  brkpt <- list()
  for(j in 1:12){ 
    ind1 = 25-(25-2*j) ; ind2 = ind1+1   #takes col 2,3 ; 3,4 ; etc (real, and lag)
    brkpt[[j]] <- breakpoints(subsfcn[,ind1] ~ subsfcn[,ind2], data = subsfcn) 
  }
  return(brkpt)
}
strchange <- lapply(subs.all, function(x) brks(x)) #breaks is the dynlm eq with breakpoints

for(i in 1:8){
  for(j in 2:13){ 
    if(j==2){
      plot(subslog[[i]][,j], ylab = expression(log[10](Fluorescence)), type="l",
           ylim = c(range(subslog[[i]][-grep("Cycle", colnames(subslog[[i]]))])), col = i)
    } 
    lines(subslog[[i]][,j])
    abline(v=strchange[[i]][[j-1]]$breakpoints, col=j-1)
    text(c(strchange[[i]][[j-1]]$breakpoints), max(subslog[[i]][-grep("Cycle", colnames(subslog[[i]]))]), 
         as.character(j), pos=2, srt=90, cex=0.65)
  }
}
return(strchange)
}

#figuring out breakpoinntss into DF or matrix
mello <- hello(subsets)




#empty df
breaks <- data.frame(
  Breaks = rep(NA, each = 96), Breaks1 = rep(NA, each = 96),
  Breaks2 = rep(NA, each = 96), Breaks3 = rep(NA, each = 96)
)



#double lapply extracts brkpts + adding NA to length
brks <- function(x){ #function for subset breakdate est.
  brkpt <- list()
  brkpt <- lapply(x, function(y) return(y$breakpoints))
  #nobrk <- lapply(x, function(y) return(length(y$breakpoints)))
  return(brkpt)
}
breaksp <- lapply(mello, function(x) brks(x))

lengthfc <- function(s, n){ #function for setting length for lists
  lapply(s, function(t) {length(t) <- n ; t})
}
n <- max(unlist(lapply(breaksp, lengths))) #max length of all lists
brkpts <- lapply(breaksp, function(x) lengthfc(x, n)) #same length for all with NA added 

empmat <- matrix(data=NA, nrow=96, ncol=5)
empmat[1:96, 1] <- matrix(cellolength[1:96], nrow=96)
for(i in 1:8){
  #ind2 = n*i ; ind1 = ind2-(n-1)
  brkpts.unlist <- lapply(brkpts, unlist)
  ind2 = length(brkpts[[i]])*i ; ind1 = ind2-(length(brkpts[[i]])-1)
  empmat[ind1:ind2, 2:(2+(n-1))] <- matrix(brkpts.unlist[[i]], byrow = T, nrow=length(brkpts[[i]]))   #[ind1:ind2]
}

empmat






#changing data to fit fcn
library(stringr)
mdata <- melt(df, id=c("Cycles"))
mdata <- mdata %>%
  dplyr::mutate(TargetName = "A") %>% 
  dplyr::mutate(cat = str_extract(variable, "[A-Z]")) %>% 
  mutate_each(funs(chartr("ABCDEFGH", "12345678", .)), cat) %>% 
  dplyr::mutate(new_idn = str_extract(variable, "\\d+$")) %>% 
  dplyr::mutate(new_id = "KW") %>% 
  dplyr::mutate(SampleID = paste0(new_id, cat, "_", new_idn)) %>% 
  select(starts_with("Cycle"), starts_with("SampleID"), starts_with("TargetName"), 
         starts_with("value"), -starts_with("new"), -starts_with("variable")) %>% 
  rename(dRn = value)
subsA  <- singtarget.list(mdata, "A") #works for single
savetarget.list(mdata) #works for general
load(file = "targ_A.Rda") #loads in as tst


######FUNCTION FIX THIS!

#full function
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/easy")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/easy", 
                pattern =  "^targ_")
orgdata = mdata #for below
library(Hmisc) #used for Lag 

brkplot <- function(orgdata, getfiles, klag, plot=FALSE){

#currently uses log10 for Fluo. add nonlogged function
  
#defining pathway and important lengths 
  files <- getfiles #pathway file
  targnames <- unique(files) #all targets
  targnames <- gsub(".Rda", "", targnames) #remove .Rda to match later
  sampnames <- mixedsort(unique(c(orgdata$SampleID))) #sorted samplenames
  #repeat all sample names for the chosen targets in files
  tmp <- rep(sampnames, length(files))
  #creating result data.frame
  targlength <- length(targnames) #length of all targets
  replength <- length(unique(orgdata$SampleID)) #length of all of the subsets
  sublength <- length(unique(gsub("_\\d+", "", sampnames))) #length of each subset: A,B,..,E

#using package strucchange to get breakpoints for each curve (lag at 2 set)
  #foreach loads in each tester, and outputs strcchange, whole output for breakpt 
  foreach(k=1:length(files)) %do% {   #load in each file
    load(file = files[[k]]) #loaded in as tst
    subs <- unlist.genparams(tst) #tester from unlistgenparams = tst loaded in
    
    #for any log10(x) where x<=0
#    try_log10 <- function(x) tryCatch({log10(x)}, error = function(e) NA)
#    subslog <- lapply(subs, function(x) x %>% 
#                        mutate_each(funs(try_log10(.)), 2:((replength/sublength)+1))) #original log10 terms 

    subslog <- lapply(subs, function(x) x %>% 
                        mutate_each(funs(log10(.)), 2:((replength/sublength)+1))) #original log10 terms 
    subslogcb <- lapply(subslog, function(x) x %>% 
                          dplyr::select(-starts_with("Cycle")) %>% #delete cycle lag
                          mutate_each(funs(Lag(., shift=klag)))) #lagged terms
#lag terms up to klag chosen
#    lagterms <- function(klag){
#      for(k in 1:klag){
#        lagforms[[k]] <- unlist(sapply("Lag(., shift=", paste0, 1:k, ")"))
#      }
#      return(lagforms[[klag]]) #returns specific klags
#    }
#    lagschosen <- lagterms(klag)
#cannot create 3columns for the klag of 1,2,3   
#    lapply(subslog, function(x) x %>% 
#             dplyr::select(-starts_with("Cycle")) %>% 
#             mutate_each(funs(lagterms(klag=3)))) #as.formula(lagschosen[[3]])

    subslogcb <- lapply(subslogcb, function(x) setNames(x, paste0(colnames(x), "_lag"))) #colnames for lagged
    subs.all <- Map(cbind, subslog, subslogcb) #mapped for column bind for log and logcb
    subs.all <- lapply(subs.all, function(x) x %>% 
                         .[mixedsort(colnames(.))] %>% #correct order except Cycles
    #dplyr::select(noquote(mixedsort(colnames(.)))) %>% #correct order except for Cycles
                         dplyr::select(starts_with("Cycle"), everything())) #correct order
    #subs.all <- lapply(subs.all, function(x) x %>% select(-starts_with("Cycle_")))
    brks <- function(subsfcn){ #function for subset breakdate est.
      brkpt <- list()
      for(j in 1:(replength/sublength)){ 
        ind1 = ((2*(replength/sublength))-(2*(replength/sublength)-2*j)) ; #25-(25-2*j)
        ind2 = ind1+1   #takes col 2,3 ; 3,4 ; etc (real, and lag) 
        brkpt[[j]] <- tryCatch({ #tryCatch any NAs for which means no breakpoints
          breakpoints(subsfcn[,ind1] ~ subsfcn[,ind2], data = subsfcn) 
      }, error=function(e) {
          return(list(breakpoints=NA))}) #return NA for breakpoints
      }
      return(brkpt)
    }
    #strchange is the whole output of breakpoint examination
    strchange <- lapply(subs.all, function(x) brks(x)) #breaks is the dynlm eq with breakpoints  
    
    #sort empmat first, so that we know how mmany breaks n to use
    #double lapply extracts brkpts + adding NA to length
    brkso <- function(x){ #function for extracting ONLY subset breakdate
      brkpt <- list()
      brkpt <- lapply(x, function(y) return(y$breakpoints))
      #nobrk <- lapply(x, function(y) return(length(y$breakpoints)))
      return(brkpt)
    }
    
    #getting only the breakpoints but with differing lengths
    breaksp <- lapply(strchange, function(x) brkso(x))
    
    #function to make lengths the same
    lengthfc <- function(s, n){ #function for setting length for lists
      lapply(s, function(t) {length(t) <- n ; t})
    }
    #nmax is length of most number of breakpoints
    nmax <- max(unlist(lapply(breaksp, lengths))) #max length of all lists
    #table of the breaks in for each subset
    nlength <- unlist(lapply(breaksp, lengths))
    #nlength above will call NA as length = 1 ; change this to NA for any NAs
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        if(lapply(breaksp, is.na)[[i]][[j]] == "TRUE"){
          nlength[[((replength/sublength)*i-((replength/sublength)-j))]] <- NA
        }
      }
    }
    brkpts <- lapply(breaksp, function(x) lengthfc(x, nmax)) #same length for all with NA added 

#create matrix out with ID and breakpoints 
    #empty matrix
    empmat <- matrix(data=NA, nrow=replength, ncol=(nmax+1))
    empmat[1:replength, 1] <- matrix(nlength[1:replength], nrow=replength) #number of breaks (1st column)
    
    #break dates for the i points into matrix, columns 2 to max in empmat
    for(i in 1:length(breaksp)){ #first col is num brks, 2 to nmax are breakpoints
      #ind2 = n*i ; ind1 = ind2-(n-1)
      brkpts.unlist <- lapply(brkpts, unlist)
      ind2 = length(brkpts[[i]])*i ; ind1 = ind2-(length(brkpts[[i]])-1)
      empmat[ind1:ind2, 2:(2+(nmax-1))] <- matrix(brkpts.unlist[[i]], byrow = T, 
                                               nrow=length(brkpts[[i]]))   #[ind1:ind2]
    }
    #empmat[1:replength, 1] <- matrix(sampnames, nrow=replength) from double to character
    
    #given break names for the 
    breaknames <- unlist(sapply("Breaks", paste0, 1:20, simplify=T)) #breakpoint names to 20 (chosen cutoff)
    empmat.df <- data.frame(empmat) %>% setnames(., c("Breaks", breaknames[1:nmax])) #names to nmax 
    #rename(Breaks = X1, Break1 = X2, Break2 = X3, Break3 = X4, Break4 = X5)
    
    #empty data frame with target info
    res0 <- data.frame(
      TargetName = rep(targnames, each = replength), SampleID = rep(sampnames, targlength),
      Group = gsub("_\\d+", "", tmp))
    breaks.mat <- cbind(res0, empmat.df)
  }
  
#plotting the strcchange break observations  
  if(plot){
    for(i in 1:sublength){
      for(j in 2:((replength/sublength)+1)){ #skip first column since cycle
       # if(j==2){ #simple connection of dots, then ablines for break obs 
          plot(subslog[[i]][,j], ylab = expression(log[10](Fluorescence)), type="l",
               ylim = c(range(subslog[[i]][-grep("Cycle", colnames(subslog[[i]]))])), col = i)
         
        #lines(subslog[[i]][,j])
        abline(v=strchange[[i]][[j-1]]$breakpoints, col=1)
        #text(c(strchange[[i]][[j-1]]$breakpoints), max(subslog[[i]][-grep("Cycle", colnames(subslog[[i]]))]), 
          #   as.character(j), pos=2, srt=90, cex=0.65)
    ind2 = length(brkpts[[i]])*i ; ind1 = ind2-(length(brkpts[[i]])-(j-1))
    legend("bottomright", as.character(breaks.mat[ind1,"SampleID"]), bty = "n", cex=0.75)
      }
    }
  }
  return(breaks.mat)
} 
  
brkplot(orgdata, getfiles, klag=5, plot=FALSE)
#larger data set
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont", 
                pattern =  "^targ_")
load(file = getfiles[[1]])
brkplot(miRcompData2, getfiles, klag=1, plot=FALSE)

####applying LSTAR model thru function
source("GAPDH.SO/lag_gen.R")
library(dynlm) ; library(Hmisc)

#add get lag formulae in there? then no source needed for lag gen

aiclag <- function(subs, log=TRUE){
#finding the lag term length that will be used by plotting AIC
fitsres <- list() ; fitsresaic <- list()

#log10 conversion
if(log){
  subslog <- lapply(subs, function(x) x %>% 
                    mutate_each(funs(log10(.)), 2:13)) #original log10 terms 
  subs = subslog
}
#time series converion of data
  subs.ts <- lapply(subs, ts) #time series of data
  
#establish get lag formulae used for later
get_lag_formulae <- function(n, nrep){
    formulae <- list()
    for(k in 1:n){
      formulae[[k]] <- paste(nrep, "~", paste(paste0("L(", nrep, ",", 1:k,")"), collapse=" + "), collapse=" ")
    }
    return(formulae)
  }

#constructing formulas
subsnames <- lapply(LETTERS[1:8], paste0, 1:12) #list of subset names
#dynlm formula in list w/ up to 15 lags
subsformulas <- lapply(subsnames, function(x) lapply(x, get_lag_formulae, n=15)) 

for(j in 1:8) fitsres[[j]] <- lapply(subsformulas[[j]], function(x) lapply(x, 
                                            function(f) dynlm(formula = as.formula(f), 
                                              data=subs.ts[[j]]))) #output results of dynlm
for(j in 1:8) fitsresaic[[j]] <- sapply(fitsres[[j]], 
                                            function(x) sapply(x,AIC)) #output results of dynlm

#plotting AIC values
for(j in 1:8){
  for(i in 1:12){
    if(i==1 & j==1){
      #cutoff 15 lags so x=1:15
      plot(x=1:15, y=fitsresaic[[j]][,i], type = "p", ylim = c(range(fitsresaic)), 
           ylab = "AIC", xlab = "Lagged Term Set")
      }
      points(x=1:15, y=fitsresaic[[j]][,i], col = j)
    } #plotting all replication sets' AIC
  }
return(fitsresaic)
}

aiclag(subs = subs, log=FALSE) #non-log10 values hard to see AIC comparison
aiclag(subs = subs) #log10 values easy to see dramatic improvement at lag=2


####STAR model function
##plotting LSTAR model to see fit for replication
lsubsets <- subsets #lsubsets different from subsets_log b/c subsets_log with ts 
for(i in 1:8){
  for(j in 2:13){
    lsubsets[[i]][[j]] <- log10(subsets[[i]][[j]])
  }
}
lsubsets <- subslog
lsubsets <- lapply(subslog, ts)

tsDyn::lstar(lsubsets[[1]][,2], m=2)

ff <- rowMeans(lsubsets$A[,2:13])
ff2 <- apply(lsubsets$A[,2:13], 1, median)

try.lstar2 <- tsDyn::lstar(fftest, m=2, d=1)
try.lstar <- tsDyn::lstar(ff, m=2, d=1) #embedding dimension=2, delay = 1 #mean

plot(x=3:40, try.lstar$fitted.values, type = "l") #ylim = c(5.1, 5.75), xlim=c(1,40))
for(j in 1:12){
  for(h in 1:8){
    points(x=1:40, y=subsets_log[[h]][[j]][,1], cex=0.45)
  }
} #plot points for 38 cycles given lag 2try.lstar2 <- tsDyn::lstar(ff2, m=2, d=1) #median replication set
#lines(x=1:38, try.lstar2$fitted.values, type = "l", ylim = c(5.1, 5.75), col = 2)

df_b5_log <- genparams(est=b5, listdf=lsubsets)
lines(x=1:40, y=b5_model(1:40, b=df_b5_log$params$b[1], c=df_b5_log$params$c[1],
                         d=df_b5_log$params$d[1], e=df_b5_log$params$e[1], 
                         f=df_b5_log$params$f[1]), col=2) #lines for b_5 model


#plotting residuals: do we need to forecast, if by looking, resid worse for lstar than log models?
plot(x=1:38, try.lstar$residuals, type = "l") #reduces branching effect?? doesnt match up

ff <- vector("list", 8)
for(i in 1:8){
  ff[[i]] <- rowMeans(lsubsets[[i]][, 2:13])
} #list of 8 for row means

try.lstar <- vector("list", 8)
for(i in 1:8){
  try.lstar[[i]] <- tsDyn::lstar(ff[[i]], m=2, d=1) #run lstar model through each rep
}

for(i in 1:8){
  if(i == 1){ #plot residuals/ note: in ln units
    plot(x = 1:38, y = try.lstar[[i]]$residuals, type = "l", col = 1, cex = 1.5, 
         ylab = "Residuals", xlab = "Time Series Cycle (t)")
  }
  lines(x = 1:38, y = try.lstar[[i]]$residuals, col = i, cex = 1.5)
}
legend("bottomright", c(LETTERS[1:8]), col=1:8, ncol = 4, lty = 1)
