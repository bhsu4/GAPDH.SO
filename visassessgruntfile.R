#get the files which overwrites tst 
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont", 
                pattern =  "^targ_")
files <- getfiles

#matrix with qc
res_b5 <- res_qc(getfiles, orgdata = miRcompData2, est = b5, b5dat2)
res_l5 <- res_qc(getfiles, orgdata = miRcompData2, est = l5, l5dat2)
res_b4 <- res_qc(getfiles, orgdata = miRcompData2, est = b4, b4dat2)
res_l4 <- res_qc(getfiles, orgdata = miRcompData2, est = l4, l4dat2)

save(res_b5, file = paste0("resb5" , ".Rda"))
save(res_l5, file = paste0("resl5" , ".Rda"))
save(res_b4, file = paste0("resb4" , ".Rda"))
save(res_l4, file = paste0("resl4" , ".Rda"))

#loading in matrices of parameter est
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
getmat <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont", 
              pattern = "res")

load(file = getmat[[1]]) ; load(file = getmat[[2]]) ; load(file = getmat[[3]]) ; load(file = getmat[[4]])

#creating same structure as qpcRb5
qpcRb5$ct


ind2 <- length(unique(orgdata$SampleID))*k  ; ind1 <- ind2-(length(unique(orgdata$SampleID))-1)


res_b4[1:40,]$TargetName.x
res_b4[1:40,]$ct
res_b4[1:40,]$Rsq

res_b4[1:40,]$SampleID.x

t(res_b4[1:40, c("SampleID.x", "ct")])
cbind(data.frame)
  
cbind(data.frame(TargetName = unique(res_b4[1:40,]$TargetName.x)), t(res_b4[1:40, c("ct")]))

#mixedsort the columns, then cbind with colnames = vector of sampleid

water <- res_b4[1:40, c("SampleID.x", "ct")]
sampids <- res_b4[1:40, c("SampleID.x")]
river <- order(match(sampids, paste0(hello, "_", 1:4)))
res_test <- t(water[river,]$ct) 

res_test2 <- cbind(data.frame(TargetName = unique(res_b4[1:40,]$TargetName.x)), t(as.vector(res_test)))
colnames(res_test2) <- c("TargetName", as.vector(sampids[river]))
