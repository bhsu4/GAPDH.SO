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

##written as function
target.list <- function(orgdata, target){
  
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
    
  grp.list <- list() ; repl <- list()
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
rm(grp.list)
targetatt <- target.list(miRcompData, target = targnames[1])




