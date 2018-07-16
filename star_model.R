#keeps first column cycle, log10 rest of fluorescence 
library(tidyverse) ; library(gtools) ; library(strucchange)

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
breaks <- lapply(subs.all, function(x) brks(x)) #breaks is the dynlm eq with breakpoints

for(i in 1:8){
  for(j in 2:13){ 
    if(j==2){
      plot(subslog[[i]][,j], ylab = expression(log[10](Fluorescence)), type="l",
           ylim = c(range(subslog[[i]][-grep("Cycle", colnames(subslog[[i]]))])), col = i)
    } 
    lines(subslog[[i]][,j])
    abline(v=breaks[[i]][[j-1]]$breakpoints, col=j-1)
    text(c(breaks[[i]][[j-1]]$breakpoints), max(subslog[[i]][-grep("Cycle", colnames(subslog[[i]]))]), 
         as.character(j), pos=2, srt=90, cex=0.65)
  }
}
return(breaks)
}

#figuring out breakpoinntss into DF or matrix
mello <- hello(subsets)

#eempty df
breaks <- data.frame(
  Breaks = rep(NA, each = 96), Breaks1 = rep(NA, each = 96),
  Breaks2 = rep(NA, each = 96), Breaks3 = rep(NA, each = 96)
)

#double lapply extracts brkpts
brks <- function(x){ #function for subset breakdate est.
  brkpt <- list()
  brkpt <- lapply(x, function(y) return(y$breakpoints)) 
  #nobrk <- lapply(x, function(y) return(length(y$breakpoints)))
  return(brkpt)
}
cello <- lapply(mello, function(x) brks(x))

cello[[1]]

df <- data.frame(matrix(unlist(cello), nrow=96, byrow=T))

df <- data.frame(matrix(unlist(cello), nrow=12, byrow=T),stringsAsFactors=FALSE)



leng <- function(x){
  lapply(x, function(y) return(length(y)))
}
cellolength <- unlist(lapply(cello, function(x) leng(x)))


empmat <- matrix(data=NA,nrow=96,ncol=5)
empmat[1:96, 1] <- matrix(cellolength[1:96], nrow=96)
for(i in 1:8){
  for(j in 1:12){
empmat[i, 2:(2+empmat[i, 1])] <- matrix(cello[[i]][[j]], nrow=1)
  }
}

for(i in 1:96){
empmat[i, 2:(2+empmat[i,1])] <- lapply(cello, unlist)[1][1:(1+empmat[i,1])]
}

empmat[1, 2:5] <- matrix(cello[[1]][[1]], nrow=1)

m <- matrix(NA, nrow=96, ncol=4)

for(i in 1:8){
  for(j in 1:12){
    m <- matrix(cello[[i]][[j]], nrow=1)
  }
}

8*1 - 7, 8*2 - 14, .... 

1,...12, 13

#changing data to fit fcn
mdata <- melt(df, id=c("Cycles"))
mdata <- mdata %>%
  dplyr::mutate(TargetName = "A") %>% 
  dplyr::mutate(cat = str_extract(variable, "[A-Z]")) %>% 
  mutate_each(funs(chartr("ABCDEFGH", "12345678", .)), "cat") %>% 
  dplyr::mutate(new_idn = str_extract(variable, "\\d+$")) %>% 
  dplyr::mutate(new_id = "KW") %>% 
  dplyr::mutate(SampleID = paste0(new_id, cat, "_", new_idn)) %>% 
  select(starts_with("Cycle"), starts_with("SampleID"), starts_with("TargetName"), 
         starts_with("value"), -starts_with("new"), -starts_with("variable")) %>% 
  rename(dRn = value)
subsA  <- singtarget.list(mdata, "A") #works for single
savetarget.list(mdata) #works for general
load(file = "targ_A.Rda") #loads in as tst


#overlaps text when attaching at brkpt
#create DF of it, return that DF, and paste multiple cycle number on graph

getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/easy", 
                pattern =  "^targ_")

files <- getfiles
targnames <- unique(files) #all targets
targnames <- gsub(".Rda", "", targnames) #remove .Rda to match later
sampnames <- mixedsort(unique(c(orgdata$SampleID)))
#repeat all sample names for the chosen targets in files
tmp <- rep(sampnames, length(files))
#creating result data.frame
targlength <- length(targnames) #length of all targets
replength <- length(unique(orgdata$SampleID)) #sum(lengths(tst)) #length of total reps for each target

res <- data.frame(
    TargetName = rep(targnames, each = replength), SampleID = rep(sampnames, targlength),
    Group = gsub("_." , "", tmp), Breaks = rep(NA, each = replength), 
    Break1 = rep(NA, each = replength), Break2 = rep(NA, each = replength),
    Break3 = rep(NA, each = replength)
)

res[1, 5:7] <- c(mello[[1]][[1]]$breakpoints)


meck <- list()
for(i in 1:8){
  meck[i] <- lapply(mello[[i]], function(x) length(x$breakpoints))
}

foreach(k=1:length(files)) %do% {    #(k in 1:length(files)){  #foreach (k = 1:length(files)) %do% {
  load(file = files[[k]])
  try <- unlist.genparams(tst)
  ind2 <- length(unique(orgdata$SampleID))*k  ; ind1 <- ind2-(length(unique(orgdata$SampleID))-1)
  res[ind1:ind2, 5:18] <- plot_sig(est, try)
  res[ind1:ind2, "FeatureSet"] <- rep(as.character(unique(lapply(tst[[1]], 
                                  function(x) unique(x$FeatureSet)))), replength)
  }


