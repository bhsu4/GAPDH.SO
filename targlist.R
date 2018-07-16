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
savetarget.list <- function(orgdata){
  
  targnames <- unique(c(orgdata$TargetName)) 
  controls <- c(names(which(table(orgdata$TargetName)>10000))) #controls
  if(length(controls) != 0){
    targnames <- targnames[!grepl(paste0(controls, collapse = "|"), targnames)]
  }
  else{}
  
#creating a group  
  splitgroup <- strsplit(orgdata[,"SampleID"], "_") #split into two parts: KW3_1 to KW3/1
  ind.keep <- seq(1,dim(orgdata)[1],1) #into sequence need to take 1st part
  splitgroup.unlist <- unlist(splitgroup) #double cancles out: all miRcompData2 (same length)
  mirc.gr <- splitgroup.unlist[seq(1,length(splitgroup.unlist), 2)] #group names -- keep odds ex: KW3 part (remove _1)
#adding group column  
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
