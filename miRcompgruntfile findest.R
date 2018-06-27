setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetssmall")
files <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetssmall", 
             pattern =  "^targ_")
for(k in 1:length(files)){
  load(file = files[[k]])
}
targnames <- unique(files)    # unique(c(miRcompData2$TargetName)) 
targnames <- gsub(".Rda", "", targnames)
sampnames <- sort(unique(c(miRcompData2$SampleID)))

load(file=files[[1]])

tmp <- rep(sampnames, 5)     #758
res <- data.frame(
  TargetName = rep(targnames, each = 40), 
  SampleID = rep(sampnames, 5), 
  Group = gsub("_." , "", tmp), 
  FeatureSet = rep(unique(tst[[1]][[1]]$FeatureSet), each = 40), 
  b = rep(NA, 5 * 40), 
  c = rep(NA, 5 * 40), 
  d = rep(NA, 5 * 40), 
  e = rep(NA, 5 * 40), 
  f = rep(NA, 5 * 40), 
  r.amp = rep(NA, 5 * 40), 
  dw.amp = rep(NA, 5 * 40), 
  p.amp = rep(NA, 5 * 40), 
  r.res = rep(NA, 5 * 40), 
  dw.res = rep(NA, 5 * 40), 
  p.res = rep(NA, 5 * 40), 
  rss = rep(NA, 5 * 40), 
  ct = rep(NA, 5 * 40), 
  eff = rep(NA, 5 * 40)
)


library(profvis)
profvis({
  for(k in 1:5){   #ength(files)){
    load(file = files[[k]])
    try <- unlist.genparams(tst)
    ind2 <- 40*k  ; ind1 <- ind2-39
    if((grepl(tst[[1]][[1]]$TargetName[[1]], res[,"TargetName"][ind1])) == "TRUE"){
      res[ind1:ind2, 5:18] <- plot_sig(b5, try)
    }
    else{print("NOPE")}
  }
})

for(k in 1:5){   #ength(files)){
  load(file = files[[k]])
  try <- unlist.genparams(tst)
  ind2 <- 40*k  ; ind1 <- ind2-39
  if((grepl(tst[[1]][[1]]$TargetName[[1]], res[,"TargetName"][ind1])) == "TRUE"){
    res[ind1:ind2, 5:18] <- plot_sig(b5, try)
  }
  else{print("NOPE")}
}






#rbind matrix