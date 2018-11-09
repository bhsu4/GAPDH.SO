#load in miRcompData
library(miRcompData)
data("miRcompData")
miRcompData2 <- miRcompData
miRcompData2 <- miRcompData2[!(miRcompData2$SampleID == "E1_" | miRcompData2$SampleID == "E2_"),] #delete E's

#changing GAPDH.SO data to fit function
library(stringr) ; library(dplyr)
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

#load in functions for LSTAR
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/")
source("GAPDH.SO/STAR_model.R")

#larger data set
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello", 
                pattern =  "^targ_")
library(foreach) ; library(Hmisc) ; library(dplyr) ; library(strucchange) ; library(gtools)
testdb <- brkplot(miRcompData2, getfiles, klag=1, plot=FALSE)

#AIC quality control for GAPDH.SO
dev.off()
set1aic <- aiclag(subs = subs, klag=15) 

#LSTAR cannot run with dplyr in library
detach("package:dplyr") ; library(dynlm) ; library(car)
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello", 
                pattern =  "^targ_")
tstlstarmat <- plot_lstar(miRcompData2, getfiles, klag=2, mdim=1, breakdb = testdb, plot=FALSE) #noplot
tstlstarmat2 <- plot_lstar(miRcompData2, getfiles, klag=1, mdim=1, breakdb = testdb, plot=TRUE) #plot
tstlstarmat3 <- plot_lstar(miRcompData2, getfiles, klag=1, mdim=2, breakdb = testdb, plot=TRUE) #plot

#specific noisy data set testing
#not nice dataset
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/noisy")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/noisy", 
                pattern =  "^targ_")
library(foreach) ; library(Hmisc) ; library(dplyr) ; library(strucchange) ; library(gtools)
testdb <- brkplot(miRcompData2, getfiles, klag=1, plot=FALSE)
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/noisy/targ_hsa-let-7c#_002405.Rda")
try <- unlist.genparams(tst)
detach("package:dplyr") ; library(dynlm) ; library(car)
tstlstarmat2 <- plot_lstar(miRcompData2, getfiles, klag=1, mdim=1, breakdb = testdb, plot=TRUE) #plot

#plotting GAPDH.SO
names(subsets$F)[[8]] <- 'F7'
library(foreach) ; library(Hmisc) ; library(dplyr) ; library(strucchange) ; library(gtools)
testdb2 <- brkplot(mdata, getfiles = NULL, subs_list=subsets, klag=2, kbreaks = 1, plot=FALSE)
detach("package:dplyr") ; library(dynlm) ; library(car)
tstlstarmat2 <- plot_lstar(mdata, getfiles=NULL, subs_list = subsets, klag=1, mdim=1, breakdb = testdb2, plot=TRUE) #plot

lstarmodnolog <- list()
for (j in 1:8) lstarmodnolog[[j]] <- lapply(subsets[[j]][,2:13], function(f) try(tsDyn::lstar(f, m=2, d=1)))

par(mfrow=c(1,1))
for(i in 12:12){
  if(i==12){
    plot(x=3:40, y=lstarmodnolog[[8]][[i]]$residuals, type = 'p', 
         col = 1, ylim = c(range(lapply(lstarmodnolog[[7]], function(x) unlist(x$residuals)))))
  }
  else{
    lines(x=3:40, y=lstarmodnolog[[8]][[i]]$residuals, col=i)
  }
}

b5resids <- lapply(df_b5$fits, resid)
b5resids_list <- vector('list', 8) ; resids.dw <- list()
for(i in 1:8){
  for(j in 1:12){
    ind2 = 40*j ; ind1 = (40*j)-(39)
    b5resids_list[[i]][[j]] <- b5resids[[i]][ind1:ind2]
    resids.dw[[12*(i-1)+j]] <- sum((b5resids_list[[i]][[j]]-Lag(b5resids_list[[i]][[j]], 1))^2, na.rm=TRUE)/sum(b5resids_list[[i]][[j]]^2)
  }
}
ml1 <- list() ; res1 <- list()
for(i in 1:8){
ml1[[i]] <- modlist(subsets[[i]], model=l5)
res1[[i]] <- getPar(ml1[[i]], type = "curve", cp = "cpD2", eff = "sliwin")
}
ct.l5 <- unlist(lapply(res1, function(x) x[1,]))

smoothScatter(x=resids.dw, y=ct.l5)
