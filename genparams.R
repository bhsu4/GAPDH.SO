genparams <- function(est, listdf){
  n <- length(listdf)    #unique(gsub("[[:digit:]+ | [:lower:] | \\.]","", colnames(df))))
  listmod <- list("l4", "l5", "b4", "b5")
  result <- lapply(listdf, function(k) pcrfit(k, fluo = 2:length(k), model = est, start = NULL,
                                              offset = 0, weights = NULL, verbose = TRUE))
  if(any(gsub("[[:alpha:]]","", result[[1]]$MODEL$name) == "5") == "TRUE") {
    for (k in 1:n){
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
#    h = which(listmod == result[[1]]$MODEL$name)
#    assign(paste0("DF", h), newtest)
  }
  if(any(gsub("[[:alpha:]]","", result[[1]]$MODEL$name) == "4") == "TRUE") {
    for (k in 1:n){
      if (k < 2) {
        params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
        test <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                             params[,"e"]))
      }
      params <- apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
      test[,k] <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                               params[,"e"]))
    }
    colnames(test) <- c("A", "B", "C", "D", "E", "F", "G", "H")
    newtest <- data.frame(t(test))
#    h = which(listmod == result[[1]]$MODEL$name)
#    assign(paste0("DF", h), newtest)
  }
  return(list(params=newtest, fits=result))
}

sub_genparams <- function(est, listdf){
  n <- length(listdf)    #unique(gsub("[[:digit:]+ | [:lower:] | \\.]","", colnames(df))))
  result = list() ; 
  for(i in 2:n){
    result[[i-1]] <- try(pcrfit(listdf, fluo=i, model = est, start = NULL,
                                offset = 0, weights = NULL, verbose = TRUE), silent=TRUE)
  }
  #res.mod <- sapply(result.tst, function(m) print(m["MODEL"])) #extract MODEL
  #mod.nam <- sapply(res.mod, function(n) print(n["name"]))     #extract names in MODEL
#  mod.nam <- unique(as.vector(unlist(sapply(result, 
#                                 function(m) sapply(m["MODEL"], 
#                                                    function(d) d["name"]))))) #vector unique names
#  mod <- unique(unlist(mod.nam)[!is.na(unlist(mod.nam))])
#  if(any(gsub("[[:alpha:]]","", mod) == "5") == "TRUE"){      #result[[1]]$MODEL$name) == "5"
   if(est$name == "l5" || est$name == "b5"){
    for (k in 1:(n-1)){
        params <- tryCatch({
          apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
          
        }, error = function(e){
          return(matrix(data=rep(NA,5), nrow=1, byrow=FALSE, 
                        dimnames=list(c(""),c("b", "c", "d", "e", "f"))))
        })
      if(k==1){
        test <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                             params[,"e"], params[,"f"]))
        }
      if(k>1){
        test[,k] <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                                 params[,"e"], params[,"f"]))   
        }
      }
      colnames(test) <- c(LETTERS[1:(n-1)])
      newtest <- data.frame(t(test))
      colnames(newtest) <- c("b", "c", "d", "e", "f")
    }
#  if(any(gsub("[[:alpha:]]","", mod) == "4") == "TRUE"){      #result[[1]]$MODEL$name) == "4"
   if(est$name == "l4" || est$name == "b4"){
    for (k in 1:(n-1)){
      params <- tryCatch({
        apply(result[[k]]$parMat[2,-1,drop=FALSE], c(1,2), as.numeric)
        
      }, error = function(e){
        return(matrix(data=rep(NA,4), nrow=1, byrow=FALSE, 
                      dimnames=list(c(""),c("b", "c", "d", "e"))))
      })
      if(k==1){
        test <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                             params[,"e"]))
      }
      if(k>1){
        test[,k] <- data.frame(c(params[,"b"], params[,"c"], params[,"d"], 
                                 params[,"e"]))   
      }
    }
    colnames(test) <- c(LETTERS[1:n-1])
    newtest <- data.frame(t(test))
    colnames(newtest) <- c("b", "c", "d", "e")
    }
  return(list(params=newtest, fits=result))
}

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
