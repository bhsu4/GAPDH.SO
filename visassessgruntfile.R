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

ind2 <- length(unique(orgdata$SampleID))*k  ; ind1 <- ind2-(length(unique(orgdata$SampleID))-1)

water <- res_b4[1:40, c("SampleID.x", "ct")]
sampids <- res_b4[1:40, c("SampleID.x")]
river <- order(match(sampids, paste0(hello, "_", 1:4)))
res_test <- t(water[river,]$ct) 
res_test2 <- cbind(data.frame(TargetName = unique(res_b4[1:40,]$TargetName.x)), t(as.vector(res_test)))
colnames(res_test2) <- c("TargetName", as.vector(sampids[river]))


#in function terms

for(k in 1:2){
ind2 <- length(unique(orgdata$SampleID))*k  ; ind1 <- ind2-(length(unique(orgdata$SampleID))-1)

#set up the df on first run through
if(k == 1){
  water <- res_b4[ind1:ind2, c("SampleID.x", "ct")]
  sampids <- res_b4[ind1:ind2, c("SampleID.x")]
  river <- order(match(sampids, paste0(mixedsort(gsub( "_.*$", "", res_b4[1:40, c("SampleID.x")])), "_", 1:4)))
  res_test <- t(water[river,]$ct) 
  res_test2 <- cbind(data.frame(TargetName = unique(res_b4[ind1:ind2,]$TargetName.x)), t(as.vector(res_test)))
  colnames(res_test2) <- c("TargetName", as.vector(sampids[river]))
}

else{
  water <- res_b4[ind1:ind2, c("SampleID.x", "ct")]
  sampids <- res_b4[ind1:ind2, c("SampleID.x")]
  river <- order(match(sampids, paste0(mixedsort(gsub( "_.*$", "", res_b4[1:40, c("SampleID.x")])), "_", 1:4)))
  res_test <- t(water[river,]$ct) 
  res_test3 <- cbind(data.frame(TargetName = unique(res_b4[ind1:ind2,]$TargetName.x)), t(as.vector(res_test)))
  colnames(res_test3) <- c("TargetName", as.vector(sampids[river]))
  res_resf <- rbind(res_test2, res_test3)
  }
}
#create a df, and matrix, then bind it together

res_ct <- matrix(NA, length(unique(res_b4$TargetName.x)), 40)
res_rsq <- matrix(NA, length(unique(res_b4$TargetName.x)), 40)

for(k in 1:length(unique(res_b4$TargetName.x))){
  ind2 <- length(unique(orgdata$SampleID))*k  ; ind1 <- ind2-(length(unique(orgdata$SampleID))-1)
  res_curr <- res_b4[ind1:ind2, c("SampleID.x", "ct", "Rsq")]
  sampids <- res_b4[ind1:ind2, c("SampleID.x")] #all unique sample IDs (KW)
  sampids_ord <- order(match(sampids, paste0(mixedsort(gsub( "_.*$", "", res_b4[ind1:ind2, c("SampleID.x")])), "_", 1:4)))
  ct_val <- t(res_curr[sampids_ord,]$ct) #correct order, mixed sorting
  rsq_val <- t(res_curr[sampids_ord,]$Rsq)
  #creating matrix for ct
  res_ct[k,] <- ct_val
  res_rsq[k,] <- rsq_val
  #matrix names
  rownames(res_ct) <- unique(res_b4$TargetName.x)
  rownames(res_rsq) <- unique(res_b4$TargetName.x)
  colnames(res_ct) <- as.vector(paste0(rep(paste0("KW", 1:10), each = 4), "_", 1:4))
  colnames(res_rsq) <- as.vector(paste0(rep(paste0("KW", 1:10), each = 4), "_", 1:4))
}
res_qpcr <- list(ct=res_ctf, qc=res_rsqf)
