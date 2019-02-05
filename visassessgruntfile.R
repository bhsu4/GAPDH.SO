#get the files which overwrites tst 
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
getfiles <- dir(path = "C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont", 
                pattern =  "^targ_")
files <- getfiles

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



