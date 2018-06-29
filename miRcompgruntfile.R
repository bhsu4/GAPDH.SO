library(miRcompData)
data("miRcompData")

table(miRcompData$SampleID)
E1 <- miRcompData[grep("E1_", miRcompData$SampleID), ]

for(k in 1:8){
  for(i in 1:12){
    lines(x=5:35, y=breakA[[1]][[k]][,i], col=k)
  } #breakpoints for replicates
}


str(let7a)
table(let7a$TargetName)
sum(let7a$Cycle == 1)
let7a_mat <- matrix(let7a$Rn, ncol = 42)
str(let7a_mat)
length(let7a$Rn)
let7a_KW5 <- matrix(let7a$Rn[grep("KW5", let7a$SampleID)], ncol= 4)
str(let7a_KW5)

tst <- miRcompData[grep("KW5", miRcompData$SampleID),]
dim(tst)
tst2 <- tst[grep("hsa-miR-21", miRcompData$TargetName),]
str(tst2)
table(tst2$TargetName)
tst2 <- tst[grep("hsa-miR-21_000397", tst$TargetName),]
tst3 <- tst[grep("hsa-miR-642_001592", tst$TargetName),]

dim(tst2)
plot(x=tst2$Cycle, y=tst$Rn)

targnames <- unique(c(miRcompData$TargetName)) 
sampnames <- unique(c(miRcompData$SampleID))

####

#42 samples with 1 replicate each -- specify target name 
testing <- miRcompData[grep(targnames[[5]], miRcompData$TargetName),]
testsamplist <- list()
for(k in 1:length(sampnames)){
  testsamplist[[k]] <- testing[grep(sampnames[[k]], testing$SampleID),]
}
plot(testsamplist[[1]]$Cycle, y=testsamplist[[1]]$dRn)


target <- function(){
  testing
}

#list by samples (42)
samplist <- list()
for(k in 1:length(sampnames)){
  samplist[[k]] <- miRcompData[grep(sampnames[[k]], miRcompData$SampleID),]
} 

for(k in 1:length(targnames)){
  if(unique(samplist[[k]]["SampleID"])[[1]] = targnames[[k]]){
    
  }
  
}

unique(samplist[[1]]$TargetName)
length(unique(samplist[[1]]$TargetName)) #758 unique targets regardless of sampleID

#assume samplist[[1]] = KW3_1
for(j in 1:1){ #42 length(unique(samplist[[1]]$SampleID))
  for(i in 1:length(unique(samplist[[1]]$TargetName))){ #758
    if(samplist[[1]]$TargetName == targnames[[i]]){
      lapply(samplist,function(x) grep(targnames[[i]],x)) 
      
      #temp <- data.frame(samplist[[1]])
      #temp <- temp[grep(targnames[[i]], miRcompData$TargetName), ]
      plot(x=temp$Cycle, y=temp$dRn)
    }
  }
}

for(i in 1:42){
  samplist[[i]] <- samplist[[i]][grep(targnames[[1]], samplist[[i]]["TargetName",])]
}


testing <- lapply(samplist, function(x) grep(targnames[[1]], x)) 



testing <- lapply(samplist, function(TargetName) TargetName[,grep(targnames[[1]],colnames(TargetName))])





if(unique(samplist[[k]]["SampleID"])[[1]] = targnames[[k]]){
  plot(x=samplist[[1]]$Cycle, y=samplist[[1]]$dRn)
}






####creating new mircompdata with group

miRcompData2 <- miRcompData
miRcompData2 <- miRcompData2[!(miRcompData2$SampleID == "E1_" | miRcompData2$SampleID == "E2_"),]

splitgroup <- strsplit(miRcompData2[,3], "_") #split into two parts: KW3_1 to KW3/1
ind.keep <- seq(1,dim(miRcompData2)[1],1) #into sequence need to take 1st part
splitgroup.unlist <- unlist(splitgroup) #double cancles out: all miRcompData2 (same length)
mirc.gr <- splitgroup.unlist[seq(1,length(splitgroup.unlist), 2)] #group names -- keep odds ex: KW3 part (remove _1)

mirc.order <- order(mirc.gr,decreasing = FALSE) 
new.mirc.data <- miRcompData[mirc.order,]
new.mirc.data <- cbind(miRcompData2[mirc.order, ], group = mirc.gr[mirc.order]) #both orders same
unique(new.mirc.data$group) #only 10 groups 

##### not I%n% function
#'%!in%' <- function(x,y)!('%in%'(x,y))
#excl.gr <- which(unlist(splitgroup)=="E1" | unlist(splitgroup)=="E2")
#ind.keep <- ind[ind %!in% excl.gr]

targnames[1] -> targ # target
# function for creating is of lists for each target
#ooefjsgj <- function(my.target){

for(h in 1:length(targnames)) { #repeats this 758 times for 10 by 4 lists
  targ  <- targnames[h]
grp.list <- list()
for(j in 1:length(unique(new.mirc.data$group))){ #list w/ groups of 10 KWs (primary list)
    grp.list[[j]] <- new.mirc.data[which(new.mirc.data$group == unique(new.mirc.data$group)[j] 
                                    & new.mirc.data$TargetName == targ), ] 
}

get.repl <- function(tst){ #function for secondary rep group 
  repl <- list()
  for(k in 1:length(unique(tst$SampleID))){
    if (dim(tst)==0) {print("no such combo")} else {
      #   print(k); print( unique(tst$SampleID)[k])
      repl[[k]] <- tst[which(tst$SampleID == unique(tst$SampleID)[k]),] #creating four reps within 10 reps
    } #tst is a list of 10
  }
  return(repl)
}

for(i in 1:10){
  tst[[i]] <- get.repl(grp.list[[i]])
}
save(tst, file = paste0("LOL_",targ, ".Rda"))

}

##written as function (single target specified)
singtarget.list <- function(orgdata, target){
  
  targnames <- unique(c(orgdata$TargetName)) 
  
  splitgroup <- strsplit(orgdata[,"SampleID"], "_") #split into two parts: KW3_1 to KW3/1
  ind.keep <- seq(1,dim(orgdata)[1],1) #into sequence need to take 1st part
  splitgroup.unlist <- unlist(splitgroup) #double cancles out: all miRcompData2 (same length)
  mirc.gr <- splitgroup.unlist[seq(1,length(splitgroup.unlist), 2)] #group names -- keep odds ex: KW3 part (remove _1)
  
  mirc.order <- order(mirc.gr,decreasing = FALSE) 
  ndata <- orgdata[mirc.order,]
  ndata <- cbind(orgdata[mirc.order, ], group = mirc.gr[mirc.order]) #both orders same
  unique(ndata$group) #only 10 groups 
  
#for(h in 1:length(targnames)) { #repeats this 758 times for 10 by 4 lists
    
  grp.list <- list() ; repl <- list() ; tst <- list()
    for(j in 1:length(unique(ndata$group))){
      grp.list[[j]] <- ndata[which(ndata$group == unique(ndata$group)[j]
                                   & ndata$TargetName == target), ]
    }

  #  for(k in 1:length(unique(grp.list$SampleID))){
  #    if(dim(grp.list)==0) {print("no such combo")} else{ #will print if no target
  #      repl[[k]] <- grp.list[which(grp.list$SampleID == unique(grp.list$SampleID)[k]), ]
  #    }
  #  }
  
  #return primary list of 10, and secondary list of 4
  for(i in 1:length(unique(ndata$group))){
    tst[[i]] <- get.repl(grp.list[[i]])
  }
  return(tst)
  #save(tst, file = paste0("targ_",targ, ".Rda"))
  
  #}
}
rm(grp.list) ; rm(tst)
targetatt <- singtarget.list(miRcompData2, target = targnames[1])


##written as function (saving all target)
savetarget.list <- function(orgdata){
  
  targnames <- unique(c(orgdata$TargetName)) 
  
  splitgroup <- strsplit(orgdata[,"SampleID"], "_") #split into two parts: KW3_1 to KW3/1
  ind.keep <- seq(1,dim(orgdata)[1],1) #into sequence need to take 1st part
  splitgroup.unlist <- unlist(splitgroup) #double cancles out: all miRcompData2 (same length)
  mirc.gr <- splitgroup.unlist[seq(1,length(splitgroup.unlist), 2)] #group names -- keep odds ex: KW3 part (remove _1)
  
  mirc.order <- order(mirc.gr,decreasing = FALSE) 
  ndata <- orgdata[mirc.order,]
  ndata <- cbind(orgdata[mirc.order, ], group = mirc.gr[mirc.order]) #both orders same
  unique(ndata$group) #only 10 groups 
  
for(h in 1:length(targnames)) { #repeats this 758 times for 10 by 4 lists
  target <- targnames[h]
  grp.list <- list() ; repl <- list() ; tst <- list()
  for(j in 1:length(unique(ndata$group))){
    grp.list[[j]] <- ndata[which(ndata$group == unique(ndata$group)[j]
                                 & ndata$TargetName == target), ]
  }
  
  #return primary list of 10, and secondary list of 4
  for(i in 1:length(unique(ndata$group))){
    tst[[i]] <- get.repl(grp.list[[i]])
  }
  save(tst, file = paste0("targ_", target, ".Rda"))
  }
}

savetarget.list(miRcompData2) #saves all the miRcompData as list of lists

#loads as tst, replaces
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-let-7c#_002405.Rda")

df_b5 <- sub_genparams(est=b5, listdf=tst)
df_b5 <- genparams(est=b5, listdf=subsets)

#unlisting list of lists for secondary list is dataframe for genparams
rm(listdf.tst) ; rm(tst.list) ; rm(j) ; rm(mincyc)
tst.list <- list() ; listdf.tst <- list()
for(i in 1:10){
  j=1 #set j for within list reps
  tst.list[[i]] <- list(tst[[i]][[j]]$dRn, tst[[i]][[j+1]]$dRn, 
                        tst[[i]][[j+2]]$dRn, tst[[i]][[j+3]]$dRn)
  mincyc <- min(max(tst[[i]][[j]]$Cycle), max(tst[[i]][[j+1]]$Cycle), 
                max(tst[[i]][[j+2]]$Cycle), max(tst[[i]][[j+3]]$Cycle))
  listdf.tst[[i]] <- as.data.frame(cbind(seq(1:mincyc), unlist(tst.list[[i]][[j]])[1:mincyc], 
                                                   unlist(tst.list[[i]][[j+1]])[1:mincyc],
                                                   unlist(tst.list[[i]][[j+2]])[1:mincyc], 
                                                   unlist(tst.list[[i]][[j+3]])[1:mincyc]))
}

#unlist genparams function
unlist.genparams <- function(test){
  tst.list <- list() ; listdf.tst <- list() ; repnames <- list()
  for(i in 1:length(test)){
    for(j in 1:length(test[[i]])){ #length = rep of secondary list
      mincyc <- min(unlist(lapply(test[[i]], function(k) max(k$Cycle)))) #min cyc (if diff)
      tst.list[[i]] <- sapply(test[[i]], function(x) x$dRn[1:mincyc]) #list of diff sampleIDs, but df secondary
      listdf.tst[[i]] <- as.data.frame(cbind(seq(1:mincyc), tst.list[[i]])) 
      repnames <- lapply(LETTERS[1:length(test)], paste0, 1:length(test[[i]])) #df colnames
      names(listdf.tst[[i]]) <- c("Cycle", repnames[[i]][1:length(test[[i]])])
    }
  }
  names(listdf.tst) <- LETTERS[1:length(test)]
  return(listdf.tst)
}
rm(tst.list) ; rm(listdf.tst) ; rm(repnames) ; rm(mincyc) ; rm(test) ; rm(i)
try <- unlist.genparams(tst)

result.try <- genparams(b5, try)
source("GAPDH.SO/plot_resid.R")
plot_resid <- function(listdf, params) {
#new -- rewritten plot resid function  
  for (i in 1:length(params$fits)){
    for(k in 1:(length(listdf[[i]])-1)){
      resids <- lapply(params$fits, resid)
      if(k == 1) plot(y=resids[[i]][1:length(listdf[[i]]$Cycle)], 
                      x=params$fits[[i]]$DATA$Cycles[1:length(listdf[[i]]$Cycle)], 
                      ylim=range(resids[[i]]), type="l", 
                      xlab="Cycle", ylab="Fluoresence Residual")
      if(k > 1){
        ind2 <- length(listdf[[i]]$Cycle)*k
        ind1 <- ind2-(length(listdf[[i]]$Cycle)-1)
        lines(y=resids[[i]][ind1:ind2], x=params$fits[[i]]$DATA$Cycles[ind1:ind2], col=k)
      }
    }
    title(main= paste(names(params$fits[i]), sub=params$fits$A$MODEL$name, sep = ", "))
  }
}

result.tryresids <- lapply(result.try$fits, resid)
plot_resid(try, result.try)



##next step: re-write function / create new one: all same gene target will be together


#within replication residuals
load("C:/Users/Benjamin Hsu/Desktop/Independent Study/GAPDH.SO/targets/targ_hsa-miR-21_000397")

sub_genparams <- function(est, listdf){
  n <- length(listdf)    #unique(gsub("[[:digit:]+ | [:lower:] | \\.]","", colnames(df))))
  result = list()
  for(i in 2:n){
    result[[i-1]] <- pcrfit(listdf, fluo=i, model = est, start = NULL,
                            offset = 0, weights = NULL, verbose = TRUE)
  }
  if(any(gsub("[[:alpha:]]","", result[[1]]$MODEL$name) == "5") == "TRUE") {
    for (k in 1:(n-1)){
      if (k < 2) {
        params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
        test <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                             params[,"e"], params[,"f"]))
      }
      params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
      test[,k] <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                               params[,"e"], params[,"f"]))
    }
    colnames(test) <- c(LETTERS[1:n])
    newtest <- data.frame(t(test))
  }
  if(any(gsub("[[:alpha:]]","", result[[1]]$MODEL$name) == "4") == "TRUE") {
    for (k in 1:(n-1)){
      if (k < 2) {
        params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
        test <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                             params[,"e"]))
      }
      params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
      test[,k] <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                               params[,"e"]))
    }
    colnames(test) <- c(LETTERS[1:n-1])
    newtest <- data.frame(t(test))
  }
  return(list(params=newtest, fits=result))
}

KW1_tstresults <- sub_genparams(l4, try[[1]])

subplot_resid <- function(listdf, params){
  
  for(k in 1:(length(listdf)-1)){
    resids <- lapply(params$fits, resid)
    
    if(k == 1) plot(y=resids[[k]][1:length(listdf$Cycle)], 
                    x=params$fits[[1]]$DATA$Cycles[1:length(listdf$Cycle)], 
                    ylim=range(unlist(resids)), type="l", 
                    xlab="Cycle", ylab="Fluorescence")
    if(k > 1){
      ind2 <- length(listdf$Cycle)*k
      ind1 <- ind2-(length(listdf$Cycle)-1)*(k)-(k-1)
      lines(y=resids[[k]][ind1:ind2], x=params$fits[[1]]$DATA$Cycles[ind1:ind2], col=k)
    }
  }
}

subplot_resid(try[[1]], KW1_tstresults)

###
files <- list.files(path="path/to/dir", pattern="*.txt", full.names=T, recursive=FALSE)
lapply(files, function(x) {
  t <- read.table(x, header=T) # load file
  # apply function
  out <- function(t)
    # write to file
    write.table(out, "path/to/output", sep="\t", quote=F, row.names=F, col.names=T)
})