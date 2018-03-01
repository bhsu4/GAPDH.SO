genparams <- function(est){
  n <- length(unique(gsub("[[:digit:]+ | [:lower:] | \\.]","", colnames(df))))
  listdf <- list(A, B, C, D, E, F, G, H)
  result <- lapply(listdf, function(k) pcrfit(k, fluo = 2:13, model = est, start = NULL,
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
}
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
}

genparams(est = l4, stringsAsFactors = TRUE) #this does not work

