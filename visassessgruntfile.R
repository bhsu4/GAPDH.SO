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
load(file = "resb4.Rda") ; load(file = "resb5.Rda") ; load(file = "resl4.Rda") ; load(file = "resl5.Rda")

#creating the lists for sigmoidal models
res_qpcrl5 <- gen_ctqc(res_l5)
res_qpcrb5 <- gen_ctqc(res_b5)
res_qpcrl4 <- gen_ctqc(res_l4)
res_qpcrb4 <- gen_ctqc(res_b4)

qualityAssessment(res_qpcr)
qualityAssessment(res_qpcr, plotType="boxplot", na.rm=TRUE, label1="R-squared")
completeFeatures(res_qpcr, qcThreshold1=0.95, label1="R-sq")
boxes <- limitOfDetection(res_qpcr, qcThreshold=0.95, plotType="boxplot")
lods <- limitOfDetection(res_qpcr, qcThreshold=0.95, plotType="scatter")
lods <- limitOfDetection(res_qpcr, qcThreshold=0.95, plotType="MAplot")
titrationResponse(res_qpcrb4, qcThreshold1=0.95)
accuracy(res_qpcrb4, qcThreshold1=0.95)
boxes <- precision(res_qpcrb4, qcThreshold1=0.95, statistic="sd")

#LSTAR ct, qc scores
setwd("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targetsmcont")
load(file = "LSTAR_lag2.Rda") ; lstar_dat <- tstlstarmat3
res_qpcrlstar <- gen_ctqc(lstar_dat) #SDM
res_qpcrlstarth <- gen_ctqc(lstar_dat) #threshold method

#quality assessment
qualityAssessment(res_qpcrb4, res_qpcrlstarth, label1 = "qpcR R-squared (b4)", label2 = "LSTAR R-squared")
qualityAssessment(res_qpcrlstarth, plotType = "boxplot", na.rm=TRUE, label1 = "R-squared")
par(mfrow = c(2,1))
qualityAssessment(res_qpcrb4, res_qpcrlstarth, plotType="boxplot", na.rm=TRUE, label1="R-squared (b4)", label2 = "R-squared (LSTAR)")
#complete features
completeFeatures(res_qpcrlstarth, qcThreshold1 = 0.95, label1 = "R-sq")
#limit detection
boxes <- limitOfDetection(res_qpcrlstarth, qcThreshold=0.95, plotType="boxplot")
lods <- limitOfDetection(res_qpcrlstarth, qcThreshold=0.95, plotType="scatter")
lods <- limitOfDetection(res_qpcrlstarth, qcThreshold=0.95, plotType="MAplot")
#titration response
titrationResponse(res_qpcrlstar, qcThreshold1=0.99)
titrationResponse(res_qpcrb4, qcThreshold1=0.99, object2=res_qpcrlstarth, qcThreshold2=0.99, label1="b4", label2="lstar")
#accuracy
accuracy(res_qpcrlstar, qcThreshold1 = 0.9)
boxes <- precision(res_qpcrlstarth, qcThreshold1=0.95, statistic="sd")
accuracy(res_qpcrb4, qcThreshold1=0.99, object2=res_qpcrlstarth, qcThreshold2=0.99, label1="b4", label2="lstar")


###test out other ct methods
load(file = getfiles[[1]])
subs <- unlist.genparams(tst)
for(i in 1:sublength){ #running LSTAR model
  for(j in 2:((replength/sublength)+1)){
    #results for LSTAR model 
    lstarres[[i]][[j-1]] <- tryCatch({
      tsDyn::lstar(subs[[i]][,j], m=mdim, d=klag)}, #d = lag found through AIC
      error=function(e) list(fitted.values=rep(NA, (cyclength[[i]]-(klag*mdim))), 
                             residuals=rep(NA, (cyclength[[i]]-(klag*mdim))),
                             model.specific=list(coefficients = rep(NA, (4+2*mdim)))))
    #if error, output NA, (no fit) w/ fitted values and residual rep length times minus klag
  }
}
lstarres[[1]][[1]]$coefficients
lstarres[[1]][[1]]$model




coef <- lstarres[[1]][[2]]$coefficients

((coef["const.L"]-coef["const.H"])+ (coef["phiL.1"]-coef["phiH.1"])* Xt_1 +
  (coef["phiL.2"]-coef["phiH.2"])* Xt_2 )*(1/(1+exp(coef["gamma"]*(Xt_1-coef["th"])))) +
coef["const.H"] + coef["phiH.1"]*Xt_1 + coef["phiH.2"]*Xt_2

Xt_1 = subs[[1]][,3][2:39]
Xt_2 = subs[[1]][,3][1:38]

((coef["const.L"]-coef["const.H"])+ (coef["phiL.1"]-coef["phiH.1"])* Xt_1 +
  (coef["phiL.2"]-coef["phiH.2"])* Xt_2 )*(1/(1+exp(coef["gamma"]*(Xt_1-coef["th"])))) +
coef["const.H"] + coef["phiH.1"]*Xt_1 + coef["phiH.2"]*Xt_2

((coef["const.L"] + coef["phiL.1"]*Xt_1 + coef["phiL.2"]*Xt_2))*(1/(1+exp(coef["gamma"]*(Xt_1-coef["th"])))) +
  ((coef["const.H"] + coef["phiH.1"]*Xt_1 + coef["phiH.2"]*Xt_2))*(1-1/(1+exp(coef["gamma"]*(Xt_1-coef["th"]))))




coef["gamma"]
coef["th"]


