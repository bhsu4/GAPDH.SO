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

}

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
  rename(Rn = value)
subsA  <- singtarget.list(mdata, "A") #works!

savetarget.list(mdata)
load(file = "targ_A.Rda")


#overlaps text when attaching at brkpt
#create DF of it, return that DF, and paste multiple cycle number on graph


text(11, max(subslog[[1]][-grep("Cycle", colnames(subslog[[1]]))]), "B", pos=2, srt=90, cex=0.65)
max(subslog[[1]][-grep("Cycle", colnames(subslog[[1]]))])


targetatt <- singtarget.list(df, target = "A1") 

subletters <- LETTERS[1:8]

#if subsets is a list vs if mircmpdata2 is df
 
replength <- data.frame(lapply(subsets, nrow))

res <- data.frame(
  #target categories
  TargetName = rep(subletters, each = replength), SampleID = rep(sampnames, targlength), 
  Group = gsub("_." , "", tmp), FeatureSet = rep(NA, each = replength), 
  #parameter est for 4 parm
  b = rep(NA, targlength * replength), c = rep(NA, targlength * replength), 
  d = rep(NA, targlength * replength), e = rep(NA, targlength * replength),
  #dw statistics
  r.amp = rep(NA, targlength * replength), dw.amp = rep(NA, targlength * replength), p.amp = rep(NA, targlength * replength), 
  r.res = rep(NA, targlength * replength), dw.res = rep(NA, targlength * replength), p.res = rep(NA, targlength * replength), 
  #rss and getPar statistics
  rss = rep(NA, targlength * replength), ct = rep(NA, targlength * replength), eff = rep(NA,targlength * replength)
)

res <- data.frame(
  Target = rep(subletters[1:8], each = yoo[1:8])
)

meh = 
data.frame(target = rep(subletters[1:8], each=meh))


length(subsets$A)



















