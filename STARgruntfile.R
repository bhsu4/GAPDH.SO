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
source("GAPDH.SO/star_model")

#larger data set
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello", 
                pattern =  "^targ_")
library(foreach) ; library(Hmisc)
testdb <- brkplot(miRcompData2, getfiles, klag=2, plot=FALSE)

#AIC quality control for GAPDH.SO
dev.off()
set1aic <- aiclag(subs = subs, klag=15) 

#LSTAR cannot run with dplyr in library
detach("package:dplyr") ; library(dynlm) ; library(car)
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/mello", 
                pattern =  "^targ_")
tstlstarmat <- plot_lstar(miRcompData2, getfiles, klag=2, testdb = testdb, plot=FALSE) #noplot
tstlstarmat2 <- plot_lstar(miRcompData2, getfiles, klag=2, testdb = testdb, plot=TRUE) #plot