brkplot <- function(orgdata, getfiles, klag, kbreaks = NULL, plot=FALSE){

#currently uses log10 for Fluo. add nonlogged function
  
#defining pathway and important lengths 
  files <- getfiles #pathway file
  targnames <- unique(files) #all targets
  targnames <- gsub(".Rda", "", targnames) #remove .Rda to match later
  sampnames <- mixedsort(unique(c(orgdata$SampleID))) #sorted samplenames
  #repeat all sample names for the chosen targets in files
  tmp <- rep(sampnames, length(files))
  #creating result data.frame
  targlength <- length(targnames) #length of all targets
  replength <- length(unique(orgdata$SampleID)) #length of all of the subsets
  sublength <- length(unique(gsub("_\\d+", "", sampnames))) #length of each subset: A,B,..,E

#using package strucchange to get breakpoints for each curve (lag at 2 set)
  #foreach loads in each tester, and outputs strcchange, whole output for breakpt 
  breaksdb <- list() #return this as data.frame
  foreach(k=1:length(files)) %do% {   #load in each file
    load(file = files[[k]]) #loaded in as tst
    subs <- unlist.genparams(tst) #tester from unlistgenparams = tst loaded in
    #subslog <- lapply(subs, function(x) x %>% 
    #                    mutate_each(funs(log10(.)), 2:((replength/sublength)+1))) #original log10 terms 
    subscb <- lapply(subs, function(x) x %>% 
                      dplyr::select(-starts_with("Cycle"))) #delete cycle lag
#lag terms up to klag chosen
#    lagterms <- function(klag){
#        lagforms <- list()
#      for(k in 1:klag){
#        lagforms[[k]] <- unlist(sapply("Lag(x, shift=", paste0, 1:k, ")"))
#      }
#      return(lagforms[[klag]]) #returns specific klags
#    }
#    lagschosen <- lagterms(klag)

#lagchosen is a list of lists with primary list = lags and secondary list = values with lag(x)
  lagschosen <- list()
  lagterms <- function(x, nlag){
    y = Lag(x, shift=nlag) #specifiy nlag
    return(y)
  }
  for(i in 1:klag){ #list of lists
    lagschosen[[i]] <- lapply(subscb, function(x) lapply(x, lagterms, nlag=i)) #returns list for up to klag
  }
  #each subsets cycle length
  cyclength <- lapply(lapply(subs, function(x) x %>% select(starts_with("Cycle"))), nrow) 
  #list of 5 lags, each list with 10 reps, and each rep is df with 4 columns 
  lagschosen <- lapply(lagschosen, function(x) lapply(x, data.frame))
  
  #lagschosenn will specify the row which is also the klag specification 
  #lagschosenn <- unlist(lapply(lagschosen, function(x) lapply(x, namechanges)))
  #namechanges <- function(x){
  # print(sum(is.na(x))/(replength/sublength))
  #}
  
##extract the data, by unlist; Map accumulated df first before adding names
  #cyclength = length of each subset ; cyclengtht = total length
  cyclengtht = lapply(cyclength, "*", (replength/sublength))
  #cumcyclength is cumulative across subsets
  cumcyclength = cumsum(cyclengtht) #accumulative 
  #create lag names
  #for(i in 1:sublength){
  repnames <- unlist(sapply(LETTERS[1:sublength], paste0, 
                      1:(replength/sublength), simplify=T)) #subset rep names
  lagsnames <- unlist(sapply("Lag", paste0, 1:klag)) #Lag1, ..., Lagk names
  replags <- sapply(repnames, paste0, lagsnames) #[1:klag][[1:klag*reps*subs]] names for klags
  #}
  #laglength is the klag*reps in each subset
  laglength = length(replags)/sublength #klags * reps in each subset
  #subslag is each subset's 1:klag unlisted
  #accsubs is accumulated subsets' lags
  subslag <- list() ; accsubs <- list() 
  for(i in 1:sublength){ #number of subsets
    for(j in 1:klag){ #klag
      #each subset w/ different length ; ind1 and ind2 capture correct intervals of each subset
      ind2 = cumcyclength[[i]] ; ind1 = ind2-((cyclengtht[[i]]*(replength/sublength))- 1)
      #matrix(unlist(lagschosen[[j]])[ind1:ind2], ncol = 4) #does one group in matrix
      #subslag is subset A's (B's,..etc) unlisted klags of data (list of klag=5)
      subslag[[j]] <- unlist(lagschosen[[j]])[grep(LETTERS[i], names(unlist(lagschosen[[j]])))] #LETTERS[i]
      #accsubs is accumulated subsets w/ df of rep's lags
      #accsubs[[i]] <- matrix(unlist(subslag), ncol = laglength)
      }
      accsubs[[i]] <- data.frame(matrix(unlist(subslag), ncol = laglength))
  }
  #sdata is sequence for names
  sdata=list()
  for(i in 1:sublength){
    for(j in 1:klag){
    #ind3, ind4 is every increment for klags*reps in subset
    ind4 = laglength*i ; ind3 = ind4-(laglength-1)  
    #ind6 = (replength/sublength)*j ; ind5 = ind6-((replength/sublength)-1)
    #change name after every subset is inserted
    #colnames(accsubs[[i]])[ind5:ind6] <- replags[ind3:ind4][seq(j, laglength, 5)]  
    sdata[[j]] <- seq(j, laglength, klag) #sdata is sequential of colnames in correct order
    colnames(accsubs[[i]]) <- replags[ind3:ind4][unlist(sdata)] #5=klag -- colnames to accumulated
    }
  }
    names(accsubs) <- LETTERS[1:sublength] #names for each dataframe within accumulated list
    #ggno <- Map(cbind, lagschosen[[1]], lagschosen[[2]]) #not automated
    #attach lags with original subs data
    subs.all <- Map(cbind, subs, accsubs) #matched cbind for subs and subs w/ lag
    #paste(paste0("lagschosen", "[[", 1:5,"]]"), collapse=", ") #lagschosen[[k]] with , ; as.formula didn't work
    #subscb <- lapply(subscb, function(x) setNames(x, paste0(colnames(x), "_lag"))) #colnames for lagged
    #subs.all <- Map(cbind, subs, subscb) #mapped for column bind for log and logcb
    subs.all <- lapply(subs.all, function(x) x %>% 
                         .[mixedsort(colnames(.))] %>% #correct order except Cycles
    #dplyr::select(noquote(mixedsort(colnames(.)))) %>% #correct order except for Cycles
                         dplyr::select(starts_with("Cycle"), everything())) #correct order
    #subs.all <- lapply(subs.all, function(x) x %>% select(-starts_with("Cycle_")))
    
 ##finding the breakpoints for data using struccchange
    brk_lag_formulae <- function(nrep){ #list length(subs) w/ list of length(subsrep)
      formulae <- list()
        formulae <- paste(nrep, "~", paste(paste0(nrep, "Lag", 1:klag), 
                                collapse=" + "), collapse=" ") #generates A1 ~ A1Lag1 + A1Lag2 + A1Lagk
      return(formulae)
    }
    subsnames <- lapply(LETTERS[1:sublength], paste0, 1:(replength/sublength)) #list of subset names
    #formulas for breakpoints with inclusion of all lag terms
    subslagformulas <- lapply(subsnames, function(x) lapply(x, brk_lag_formulae)) 
    #brksres are breakpoint results: list of subs w/ reps inside [[]][[]]$breakpoints
    brksres <- list() #brksres is the whole output of breakpoint examination
    for(j in 1:sublength) brksres[[j]] <- lapply(subslagformulas[[j]], 
                                            function(x) breakpoints(formula = as.formula(x), 
                                              breaks = kbreaks, data=subs.all[[j]])) #output for breakpoints
  #**ERROR may show: must have minimum segment size > number of regressors (make klag smaller)

    #sort empmat first, so that we know how mmany breaks n to use
    #double lapply extracts brkpts + adding NA to length
    brkso <- function(x){ #function for extracting ONLY subset breakdate
      brkpt <- list()
      brkpt <- lapply(x, function(y) return(y$breakpoints))
      #nobrk <- lapply(x, function(y) return(length(y$breakpoints)))
      return(brkpt)
    }
    
    #getting only the breakpoints but with differing lengths
    breaksp <- lapply(brksres, function(x) brkso(x))
    
    #function to make lengths the same
    lengthfc <- function(s, n){ #function for setting length for lists
      lapply(s, function(t) {length(t) <- n ; t})
    }
    #nmax is length of most number of breakpoints
    nmax <- max(unlist(lapply(breaksp, lengths))) #max length of all lists
    #table of the breaks in for each subset
    nlength <- unlist(lapply(breaksp, lengths))
    #nlength above will call NA as length = 1 ; change this to NA for any NAs
    for(i in 1:sublength){
      for(j in 1:(replength/sublength)){
        if(lapply(breaksp, is.na)[[i]][[j]] == "TRUE"){
          nlength[[((replength/sublength)*i-((replength/sublength)-j))]] <- NA
        }
      }
    }
    brkpts <- lapply(breaksp, function(x) lengthfc(x, nmax)) #same length for all with NA added 

  ##create matrix out with ID and breakpoints 
    #empty matrix
    empmat <- matrix(data=NA, nrow=replength, ncol=(nmax+1))
    empmat[1:replength, 1] <- matrix(nlength[1:replength], nrow=replength) #number of breaks (1st column)
    
    #break dates for the i points into matrix, columns 2 to max in empmat
    for(i in 1:length(breaksp)){ #first col is num brks, 2 to nmax are breakpoints
      #ind2 = n*i ; ind1 = ind2-(n-1)
      brkpts.unlist <- lapply(brkpts, unlist)
      ind2 = length(brkpts[[i]])*i ; ind1 = ind2-(length(brkpts[[i]])-1)
      empmat[ind1:ind2, 2:(2+(nmax-1))] <- matrix(brkpts.unlist[[i]], byrow = T, 
                                               nrow=length(brkpts[[i]]))   #[ind1:ind2]
    }
    #empmat[1:replength, 1] <- matrix(sampnames, nrow=replength) from double to character
    
    #given break names for the 
    breaknames <- unlist(sapply("Breaks", paste0, 1:20, simplify=T)) #breakpoint names to 20 (chosen cutoff)
    empmat.df <- data.frame(empmat) %>% setnames(., c("Breaks", breaknames[1:nmax])) #names to nmax 
    #rename(Breaks = X1, Break1 = X2, Break2 = X3, Break3 = X4, Break4 = X5)
    
    #empty data frame with target info
    res0 <- data.frame(
      TargetName = rep(targnames, each = replength), SampleID = rep(sampnames, targlength),
      Group = gsub("_\\d+", "", tmp))
    res0$TargetID = seq(1, replength*length(files), 1) #each subset's replication is unique ID
    #breaks 1st normalization with repetitive columns
    breaks.mat <- cbind(res0, empmat.df)
    
  ##database for breakpoints w/ foreign key to TargetID
    #non-NA breakpoints extracted by targetname 
    nonares0 <- res0[which(res0==targnames[[k]]),][which(is.na(nlength) == FALSE),]
    #targID is targetID, breakseq is 1:n for number of breaks
    targID <- list() ; breakseq <- list()
      for(i in 1:nrow(nonares0)){ #sum of the number of targetIDs to input
        #targID is targetID repeated n times equal to number of breakpoints
        targID[[i]] <- rep(nonares0$TargetID[i], nlength[which(is.na(nlength)==FALSE)][i])
        breakseq[[i]] <- 1:length(targID[[i]]) #if 2 breaks, 1:2; 3 breaks, 1:3
      }
      breakspun <- unlist(breaksp)[which(is.na(unlist(breaksp))==FALSE)] #non-NA breakpoints
      breaksdb[[k]] <- cbind(matrix(unlist(targID), ncol=1), 
                         matrix(unlist(breakseq), ncol=1), 
                         matrix(breakspun, ncol=1)) #database for breakpoints
  }
  finaldb <- data.frame(do.call(rbind, breaksdb))
  colnames(finaldb) <- c("TargetID", "BreakSeq", "Breakpoints")
  finals <- list(Targets = res0, Breaks = finaldb)
  
#plotting the strcchange break observations  
  if(plot){
    for(i in 1:sublength){
      for(j in 2:((replength/sublength)+1)){ #skip first column since cycle
       # if(j==2){ #simple connection of dots, then ablines for break obs 
          plot(subs[[i]][,j], ylab = "Fluorescence", type="l",
               #ylim = c(range(subs[[i]][-grep("Cycle", colnames(subs[[i]]))])), col = i)
                ylim = c(range(subs[[i]][[j]])), col=i)
        #lines(subslog[[i]][,j])
        abline(v=brksres[[i]][[j-1]]$breakpoints, col=1)
        #text(c(strchange[[i]][[j-1]]$breakpoints), max(subslog[[i]][-grep("Cycle", colnames(subslog[[i]]))]), 
          #   as.character(j), pos=2, srt=90, cex=0.65)
    ind2 = length(brkpts[[i]])*i ; ind1 = ind2-(length(brkpts[[i]])-(j-1))
    legend("bottomright", as.character(breaks.mat[ind1,"SampleID"]), bty = "n", cex=0.75)
      }
    }
  }
  return(finals)
} 

aiclag <- function(subs, klag){
  #finding the lag term length that will be used by plotting AIC
  fitsres <- list() ; fitsresaic <- list()

  #defining important lengths
  targnames <- unique(files) #all targets
  targnames <- gsub(".Rda", "", targnames) #remove .Rda to match later
  sampnames <- mixedsort(unique(c(orgdata$SampleID))) #sorted samplenames
  #repeat all sample names for the chosen targets in files
  tmp <- rep(sampnames, length(files))
  #creating result data.frame
  targlength <- length(targnames) #length of all targets
  replength <- length(unique(orgdata$SampleID)) #length of all of the subsets
  sublength <- length(unique(gsub("_\\d+", "", sampnames))) #length of each subset: A,B,..,E

  #time series converion of data
    subs.ts <- lapply(subs, ts) #time series of data
  
  #establish get lag formulae used for later
  get_lag_formulae <- function(n, nrep){
      formulae <- list()
      for(k in 1:n){
        formulae[[k]] <- paste(nrep, "~", paste(paste0("L(", nrep, ",", 1:k,")"), 
                         collapse=" + "), collapse=" ") #generates L(x,1) + L(x,2) + ..L(x,k) 
      }
      return(formulae)
    }

  #constructing formulas
  subsnames <- lapply(LETTERS[1:sublength], paste0, 1:(replength/sublength)) #list of subset names
  #dynlm formula in list w/ up to 15 lags
  subsformulas <- lapply(subsnames, function(x) lapply(x, get_lag_formulae, n=klag)) 
  #dynlm formulas and AIC values
  for(j in 1:sublength) fitsres[[j]] <- lapply(subsformulas[[j]], function(x) lapply(x, 
                                            function(f) dynlm(formula = as.formula(f), 
                                              data=subs.ts[[j]]))) #output results of dynlm
  for(j in 1:sublength) fitsresaic[[j]] <- sapply(fitsres[[j]], 
                                            function(x) sapply(x,AIC)) #output results of dynlm
  #plotting AIC values (all)
  for(j in 1:sublength){
    for(i in 1:(replength/sublength)){
      if(i==1 & j==1){
        #cutoff 15 lags so x=1:15
        plot(x=1:klag, y=fitsresaic[[j]][,i], type = "p", ylim = c(range(fitsresaic)), 
             ylab = "AIC", xlab = "Lagged Term Set")
      }
        points(x=1:klag, y=fitsresaic[[j]][,i], col = j)
      } #plotting all replication sets' AIC
  }
  #individual subsets AIC plots
  for(i in 1:sublength){
    for(j in 1:(replength/sublength)){
      if(j==1){
        plot(x=1:klag, y=fitsresaic[[i]][,j], type="p", ylim = c(range(fitsresaic[[i]])),
             ylab = "AIC", xlab = "Lagged Term Set")
      }
        points(x=1:klag, y=fitsresaic[[i]][,j], col=j)
        splineres <- smooth.spline(rep(1:klag, (replength/sublength)), 
                                   y=as.vector(fitsresaic[[i]]), spar=0.35)
        lines(splineres)
    }
    legend("topright", c(subsnames[[i]]), col=1:4, lty=1, cex=0.75)
  }
  return(fitsresaic)
}

plot_lstar <- function(orgdata, getfiles, klag, mdim, breakdb, plot=FALSE){
##getting the LSTAR model parameters and coefficients
  #targets
  files <- getfiles
  #defining important lengths
  targnames <- unique(files) #all targets
  targnames <- gsub(".Rda", "", targnames) #remove .Rda to match later
  sampnames <- mixedsort(unique(c(orgdata$SampleID))) #sorted samplenames
  #repeat all sample names for the chosen targets in files
  tmp <- rep(sampnames, length(files))
  #creating result data.frame
  targlength <- length(targnames) #length of all targets
  replength <- length(unique(orgdata$SampleID)) #length of all of the subsets
  sublength <- length(unique(gsub("_\\d+", "", sampnames))) #length of each subset: A,B,..,E
  
 ###Setting Up LSTAR Output Matrix  
  res <- data.frame(
    #target categories
    TargetName = rep(targnames, each = replength), SampleID = rep(sampnames, targlength), 
    Group = gsub("_." , "", tmp), FeatureSet = rep(NA, each = replength), 
    #parameter est for LSTAR
    #const.L = rep(NA, targlength * replength), phi.L1 = rep(NA, targlength * replength), 
    #const.H = rep(NA, targlength * replength), phi.H1 = rep(NA, targlength * replength),
    #gamma = rep(NA, targlength * replength), th = rep(NA, targlength * replength),
    #dw statistics
    r.amp = rep(NA, targlength * replength), 
    dw.amp = rep(NA, targlength * replength), 
    p.amp = rep(NA, targlength * replength), 
    r.res = rep(NA, targlength * replength), 
    dw.res = rep(NA, targlength * replength), 
    p.res = rep(NA, targlength * replength), 
    #rss and getPar statistics
    rss = rep(NA, targlength * replength), 
    rssgrey = rep(NA, targlength * replength), 
    ct = rep(NA, targlength * replength)
  )

for(k in 1:length(files)){
  #loading in data sets to get subs
  load(file = files[[k]]) #loaded in as tst
  subs <- unlist.genparams(tst) #tester from unlistgenparams = tst loaded in  
  cyclength <- unlist(lapply(subs, nrow))
  #running lstar function w/ certain lag
  lstarres <- vector("list", sublength) #empty list of subs: ABCD,..etc
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
    #res1 is the list of subsets' fitted values and residuals
  res1 <- list()
  #matrix with fitted values and residuals as two columned vector
  for(i in 1:sublength){ #list of dataframes for each subset
  mat <- matrix( NA, nrow = cyclength[[i]], ncol = (2*(replength/sublength)+1))
  mat[,1] <- 1:cyclength[[i]] #different cycle lengths
    for(j in 1:(replength/sublength)){
      vector <- c(lstarres[[i]][[j]]$fitted.values, lstarres[[i]][[j]]$residuals)
      ind1 = j*2 ; ind2 = ind1+1 #gets columns 2,3 ; 3,4; 5,6
      mat[(1+(klag*mdim)):cyclength[[i]], ind1:ind2] <- matrix(vector, ncol=2) #fitted and residuals
      }
  res1[[i]] <- data.frame(mat) #df in dataframe
  }

##plotting the LSTAR model w/ actual points
 #two graphs on top each other
  par(mfrow=c(2,1))
  par(oma=c(4,4,4,4),mar=c(0.25,0.25,0,0))
#start of plot
  if(plot){
  ##REF PART 2: Clarifying Specific targetdb and breakdb
    #specific target chosen from k files
    spectarg <- breakdb$Targets[which(breakdb$Targets$TargetName == targnames[[k]]),]
    #specific break points for target
    specbreak <- breakdb$Breaks[breakdb$Breaks$TargetID %in% spectarg$TargetID, ]
  ##REF (REF PART 3): Empty Lists for Second Derivative Maximum Estimate  
    diff_df <- vector("list", 10) ; diff2_df <- vector("list", 10) 
    
  ###plotting fitted curves and actual values
  for(i in 1:sublength){ #plotting the LSTAR fitted curve
    for(j in 1:(replength/sublength)){ #i for subs, j for reps
  
  ##REF PART 3: Creating First and Second Derivatives for Maximum Estimation          
    #nonsubtractable first row    
    #lstarres[[i]][[j]]$fitted.values[-1,] 
    #nonsubtractable last row
    #lstarres[[i]][[j]]$fitted.values[-nrow(lstarres[[i]][[j]]$fitted.values),]
    #NO LSTAR Model Fit      
    if(sum(is.na(lstarres[[i]][[j]]$fitted.values)) > 0 ){
      diff_df[[i]][[j]] <-  rep(NA, length(lstarres[[i]][[j]]$fitted.values)-1)
      diff2_df[[i]][[j]] <- rep(NA, length(diff_df[[i]][[j]])-1)
    }
    #LSTAR Model
    else{
      diff_df[[i]][[j]] <-  lstarres[[i]][[j]]$fitted.values[-1,]  - lstarres[[i]][[j]]$fitted.values[-nrow(lstarres[[i]][[j]]$fitted.values),]
      }
    }
  }
  ###Removing all Neg to Pos CT values
    #unlist the first derivative slopes
    unl.diff_df <- lapply(diff_df, function(x) unlist(x, recursive=TRUE))
    for(i in 1:sublength){
        for(y in 1:length(unl.diff_df[[i]])){
      if(unl.diff_df[[i]][[y]] < 0 & is.na(unl.diff_df[[i]][[y]]) == FALSE){
          unl.diff_df[[i]][[y]] <- NA #replace all neg with NA
            }
        else{}
        }
      }
    #diff_dfl is the list of diff_df with NAs
    #get output diff2_df
    diff_dfl <- vector("list", sublength)
      for(i in 1:sublength){
        for(j in 1:(replength/sublength)){
          ind2 = length(diff_df[[i]][[j]])*j ; ind1 = ind2 - (length(diff_df[[i]][[j]])-1)
          diff_dfl[[i]][[j]] <- unl.diff_df[[i]][ind1:ind2]
          diff2_df[[i]][[j]] <- diff_dfl[[i]][[j]][-1] - diff_dfl[[i]][[j]][-length(diff_dfl[[i]][[j]])]
        }
      }
  #diff2_df[[i]][[j]] <- diff_df[[i]][[j]][-1] - diff_df[[i]][[j]][-length(diff_df[[i]][[j]])]
  #list <- unlist(lstarres, recursive = FALSE) #list of 40 equation output
  #listfittedvalues <- lapply(list, function(x) x$fitted.values) #list of 40 fitted values
  #listdftst <- do.call("cbind", listfittedvalues) #different length, cannot cbind unless NA added
  diff2_dfmax <- lapply(diff2_df, function(x) lapply(x, which.max)) #cycle at max (w/o klag)
  CTmid <- function(x) (x+(klag*mdim+1)+(x+(klag*mdim+1)+1))/2 #average of cycles
  diff2_dfCT <- lapply(diff2_dfmax, function(x) lapply(x, CTmid)) #midpoint of cycle at 2nd dermax
  
  CTrepNA <- function(x){
    if(length(x) == 0){ x <- NA }
    else{x}
  }
  diff2_dfCT <- lapply(diff2_dfCT, function(x) lapply(x, CTrepNA))
  
  #End of RSS, RSS Grey
 ###ADDED THRESHOLD CT VALUE: Beginning of RSSthres, RSS Red
 ###CT with Threshold STAR
  lstarcoef <- vector("list", sublength) #empty list of subs: ABCD,..etc
  cycCT.thover <- vector("list", sublength)
  cycCT.unl <- list() 
  #lstarcoef is list of lists of all coefficients 
  #cycCT.thover is list of lists of all cycles' LSTAR fitted > LSTAR threshold
  #cycCT.unl temporary unlisted CT.thover, because list of lists didn't work with which statement
  lstarres.unl = unlist(lstarres, recursive = FALSE)
  for(i in 1:10){
    for(j in 1:4){
      indicator <- j+4*(i-1)
      lstarcoef[[indicator]] <- lstarres.unl[[indicator]]$model.specific$coefficients #coeffs
      lstarcoef[[indicator]]["th"] <- abs(lstarcoef[[indicator]]["th"])
      if(sum(is.na(lstarcoef[[indicator]])) > 0){cycCT.unl[[indicator]] <- NA}
      else{
        #cycles greater than threshold in lagged format
        cycCT.unl[[indicator]] <- which(lstarres.unl[[indicator]]$fitted.values > lstarcoef[[indicator]]["th"])
      }
    }
    ind2 = 4*i ; ind1 = ind2-(4-1)
    cycCT.thover[[i]] <- cycCT.unl[ind1:ind2] #back to list of lists
  }
  overth.diff <- function(x){ #diff takes lagged difference, rle shows lengths of lagged diff
    #ex. x = 4,5,6,7,10,13,25,26,27,28,29,30
    #lengths are 3,2,1,5 these are how many of those differences (max this)
    #values are 1,3,12,1 these are diffrences (this needs to be = 1)
    if(sum(is.na(x)) == 0){ rle(diff(x)) }
    else{ list(lengths=NA) }
  }      
  overth.diffcyc <- function(x){ #fitting conditions above
    #max lengths and choose values of 1
    #logic here is that: in cases where more fluctuations, more noise, multiple times when
    #cycles cross threshold. take the cycle that has most cycles after threshold that is 
    #higher than the threshold continuously
    if(sum(is.na(x)) == 0){ which(x$lengths == max(x$lengths) & x$values == 1) }
    else{ NA }
  }
  #difftimes are rle (run length encoding) lagged diff
  difftimes <- lapply(cycCT.thover, function(x) lapply(x, overth.diff))
  #cycle index for cycle in list of those > threshold
  index.length <- lapply(difftimes, function(x) lapply(x, overth.diffcyc))
  index.unl <- unlist(index.length, recursive = FALSE)
  index.unl <- lapply(index.unl, max) #if two same length of lagged diff, choose latter cycles
  #final CT list
  cycCT.th <- vector("list", sublength)
  #find the intermediate CT value (intermediate means 0.5) surpassing threshold
  for(i in 1:sublength){
    for(j in 1:(replength/sublength)){
      ind = j+(replength/sublength)*(i-1)
      if(sum(is.na(index.unl[[ind]])) > 0){ #any NAs get NA
        cycCT.th[[i]][[j]] <- NA
      }
      else if(index.unl[[ind]] > 1){ #more than one string of consec > threshold
        #ex. 1,2,3,4,5,6,7,8,9,10,.....,24 (20th term), 28,29,30,....
        #difftimes lengths: 11, 1, 5, 1, 1, 1, 11
        #          values:  1, 3, 1, 3, 1, 4, 1
        #index.unl = 7, take sum of 1:6 of the lengths = 20, so 21st is our chosen pt
        tmd <- sum(difftimes[[i]][[j]]$lengths[1:(index.unl[[ind]]-1)]) #sum of before
        cycCT.th[[i]][[j]] <- (cycCT.thover[[i]][[j]][[tmd+1]] + (cycCT.thover[[i]][[j]][[tmd+1]]-1))/2 + (mdim*klag) #actual pt after add
      }
      else{ 
        #at index = 1, just take that first point, and point-1 , find average
        cycCT.th[[i]][[j]] <- (cycCT.thover[[i]][[j]][[1]] + (cycCT.thover[[i]][[j]][[1]]-1))/2 + (mdim*klag) #actual pt after addition of mdim*klag
      }
    }
  } 
  
  ###plotting fitted curves and actual values
  for(i in 1:sublength){ #plotting the LSTAR fitted curve
    for(j in 1:(replength/sublength)){ #i for subs, j for reps
#####All NO LSTAR Models#####        
  ###PART 1: NO LSTAR Model Fits with Breakpoints as Squares
    #NO LSTAR model list of NAs           
    if(sum(is.na(lstarres[[i]][[j]]$fitted.values >0))){
      plot(x=subs[[i]][[1]], y=subs[[i]][[j+1]], type="p", 
           ylab = "Fluorescence", xlab = "Cycle", xaxt = "n")
      legend("topleft", c(colnames(subs[[i]])[[j+1]]), lty=1, cex=0.65) #legend
  ###PART 2: Breakpoints as Squares for NON-LSTAR Models
    #indij for indicator with i,j used
     indij = 40*(k-1)+(4*(i-1))+j
    #sum of the TRUE > 0 then we plot square break points
    if(sum(specbreak$TargetID %in% indij) > 0){ 
    #specific breakpoint cycle number
    specbreakat <- specbreak$Breakpoints[which(specbreak$TargetID == indij)]
    #how many breakpoints for that subset
    specbreaknlength <- length(specbreakat)
   ##All Squares: No LSTAR fit, but Breakpoints present
    for(z in 1:specbreaknlength){
     points(x=specbreak$Breakpoints[which(specbreak$TargetID == indij)][[z]], 
            y=subs[[i]][[j+1]][specbreakat[[z]]],
            cex=1, pch = 15) #actual points
        }
      }
    else{}
  ###PART 3: Estimating CT Values w/ no LSTAR Model: NO CT
    #estimated CT value
    points(x=diff2_dfCT[[i]][[j]], 
           y=lstarres[[i]][[j]]$fitted.values[diff2_dfCT[[i]][[j]]-(klag*mdim)],
           pch=17, cex=1)
  ###PART 4: Empty Residuals Since NO LSTAR
    plot(1, type="n" , xlab="n", ylab="Residuals", 
         xlim=c(0, max(subs[[i]][[1]])), ylim=c(0,1), cex.lab=0.75)        
    }
        
#####All LSTAR Models#####        
  ###PART 1: LSTAR Model Fits with Breakpoints as Squares
    #LSTAR model fits               
    else{
      plot(x=(1+(klag*mdim)):length(subs[[i]][[1]]), y=lstarres[[i]][[j]]$fitted.values, 
           ylab="Fluorescence", xlab = "", type="l", xaxt = "n",
           ylim=c(min(unlist(subs[[i]][[j+1]], lstarres[[i]][[j]]$fitted.values)),
                  max(unlist(subs[[i]][[j+1]], lstarres[[i]][[j]]$fitted.values))))
      points(x=subs[[i]][[1]], y=subs[[i]][[j+1]]) #actual points
      legend("topleft", c(colnames(subs[[i]])[[j+1]]), lty=1, cex=0.65) #legend
  ###PART 2: Breakpoints as Squares for LSTAR Models
    #indij for indicator with i,j used
     indij = 40*(k-1)+(4*(i-1))+j
    #sum of the TRUE > 0 then we plot square break points
    if(sum(specbreak$TargetID %in% indij) > 0){ 
    #specific breakpoint cycle number
    specbreakat <- specbreak$Breakpoints[which(specbreak$TargetID == indij)]
    #how many breakpoints for that subset
    specbreaknlength <- length(specbreakat)
    ##All Squares: LSTAR fit with n Breakpoints present
    for(z in 1:specbreaknlength){
      points(x=specbreak$Breakpoints[which(specbreak$TargetID == indij)][[z]], 
             y=lstarres[[i]][[j]]$fitted.values[(specbreakat[[z]]-klag*mdim), ], 
             cex=1, pch = 15) #actual points
        }
      }
    else{}
  ###PART 3: Estimating CT Values w/ LSTAR Model by Slope Midpoint
    ##A. SECOND DERIVATIVE METHOD: estimated CT value            
    #*Shouldn't Display Error: klag=3, Encounter error if cycle is less than 4 
    #since starts at cyc 4, no need to filter through cycles < 4 when creating breakdb
     leftbound <- lstarres[[i]][[j]]$fitted.values[diff2_dfCT[[i]][[j]]-0.5-(klag*mdim)]
     rightbound <- lstarres[[i]][[j]]$fitted.values[diff2_dfCT[[i]][[j]]+0.5-(klag*mdim)]
     avgCTfittedvalue = (leftbound+rightbound)/2
    #estimated CT value
     points(x=diff2_dfCT[[i]][[j]], y=avgCTfittedvalue, pch=17, cex=1)
    #grey box +/- 2 cycles about estimated CT 
     polygon(x = c(diff2_dfCT[[i]][[j]]-2, diff2_dfCT[[i]][[j]]+2, diff2_dfCT[[i]][[j]]+2, 
                   diff2_dfCT[[i]][[j]]-2, diff2_dfCT[[i]][[j]]-2), 
             y = c(min(par("usr")), min(par("usr")), 
                   max(par("usr")), max(par("usr")), min(par("usr"))),
             col= rgb(0,0,0,alpha=0.15))
     ##B. LSTAR THRESHOLD METHOD: estimated CT value            
     leftbound <- lstarres[[i]][[j]]$fitted.values[cycCT.th[[i]][[j]]-0.5-(klag*mdim)]
     rightbound <- lstarres[[i]][[j]]$fitted.values[cycCT.th[[i]][[j]]+0.5-(klag*mdim)]
     avgCTfittedvalue = (leftbound+rightbound)/2
     #estimated CT value
     if(cycCT.th[[i]][[j]]>(0.5+klag*mdim)){
      #if explanation: at d*m = 1, 2:40 -> earliest cycle cross threshold is 1st cycle
      #which has intermediate at 0.5. Actual cycle is 0.5 + d*m (whch is 1) = 1.5.
      #Thus, anything less than 1.5 + 2 cycle range = 3.5 cannot have RSSred
     points(x=cycCT.th[[i]][[j]], y=avgCTfittedvalue, pch=17, cex=1, col="red")
     #grey box +/- 2 cycles about estimated CT 
     polygon(x = c(cycCT.th[[i]][[j]]-2, cycCT.th[[i]][[j]]+2, cycCT.th[[i]][[j]]+2, 
                   cycCT.th[[i]][[j]]-2, cycCT.th[[i]][[j]]-2), 
             y = c(min(par("usr")), min(par("usr")), 
                   max(par("usr")), max(par("usr")), min(par("usr"))),
             col= rgb(1,0,0,alpha=0.15))}
     else if(cycCT.th[[i]][[j]]<(0.5+klag*mdim)){ }
     
  ###PART 4: Residuals Plotting
   ##A. SECOND DERIVATIVE METHOD: estimated CT value            
    #residuals
     plot(x=(1+(klag*mdim)):length(subs[[i]][[1]]), y=lstarres[[i]][[j]]$residuals,
          ylab="Residuals", xlab = "Cycle", type="p", cex.lab=0.75)
          abline(h=0)
    #estimated CT value
     points(x=diff2_dfCT[[i]][[j]], y=0, pch=17, cex=1)
    #grey box +/- 2 cycles about estimated CT 
     polygon(x = c(diff2_dfCT[[i]][[j]]-2, diff2_dfCT[[i]][[j]]+2, diff2_dfCT[[i]][[j]]+2, 
                   diff2_dfCT[[i]][[j]]-2, diff2_dfCT[[i]][[j]]-2), 
             y = c(min(par("usr")), min(par("usr")), 
                   max(par("usr")), max(par("usr")), min(par("usr"))),
             col= rgb(0,0,0,alpha=0.15))
    ##B. LSTAR THRESHOLD METHOD: estimated CT value            
     #estimated CT value
     if(cycCT.th[[i]][[j]]>(0.5+klag*mdim)){ 
     #if explanation: at d*m = 1, 2:40 -> earliest cycle cross threshold is 1st cycle
     #which has intermediate at 0.5. Actual cycle is 0.5 + d*m (whch is 1) = 1.5.
     #Thus, anything less than 1.5 + 2 cycle range = 3.5 cannot have RSSred
     points(x=cycCT.th[[i]][[j]], y=0, pch=17, cex=1, col="red")
     #grey box +/- 2 cycles about estimated CT 
     polygon(x = c(cycCT.th[[i]][[j]]-2, cycCT.th[[i]][[j]]+2, cycCT.th[[i]][[j]]+2, 
                   cycCT.th[[i]][[j]]-2, cycCT.th[[i]][[j]]-2), 
             y = c(min(par("usr")), min(par("usr")), 
                   max(par("usr")), max(par("usr")), min(par("usr"))),
             col= rgb(1,0,0,alpha=0.15))
     }
          } #else
        } #i
      } #j
    } #if(plot)

####Creating Matrix Output w/ TargID, AC, RSS, RSSGrey
 ###PART 1: Creating the LSTAR Model's RSS and RSSGrey
   #list of lists of RSS
   rsslstar <- lapply(lstarres, function(x) lapply(x, function(x) sum(x$residuals^2)))
   #vector of 1:(replength) for RSS
   unlist.rsslstar <- unlist(unlist(rsslstar, recursive = FALSE)) #numeric vector
   #cycles for CT SDM
   unlistcycCT <- unlist(unlist(diff2_dfCT, recursive=FALSE))
   #cycles for CT LSTAR THRESHOLD
   unlistcycCTthr <- unlist(unlist(cycCT.th, recursive=FALSE))
   #*Possible Error: UnlistcycCT suggests 3:46, if <4 then klag-adjusted is from 0:3 (shown below)
   #ex. at 3.5, then 3.5-(klag=2) = 1.5 +/- 2 cycles (grey), which is 0:3.  
   rsslstarg <- vector('list', sublength) 
   for(i in 1:(sublength)){
     for(j in 1:(replength/sublength)){
       #combined indicator with i and j
       indij = ((replength/sublength)*(i-1))+j  
       #klag-adjusted upper cycle b/c unlistcycCT[[indij]] is actual cycles for 3:46 (klag=2)
       #unlistcycCT[[indij]] - klag to fit lstarres$residuals 1:44 formmat
       if(sum(is.na(unlistcycCT[[indij]])) > 0){
         print(paste(k, "-", indij, "no CT"))
         rsslstarg[[i]][[j]] <- NA
       }
       else if(sum(isTRUE(unlistcycCT[[indij]] > cyclength[[i]]-(klag*mdim+1))) > 0){ #strict upper bound
         uppercyc = cyclength[[i]]-(klag*mdim)
         lowercyc = ceiling(unlistcycCT[[indij]]-(klag*mdim)-2)}
       else{ 
         uppercyc = floor(unlistcycCT[[indij]]-(klag*mdim)+2)  #+2 cycles
         #klag-adjusted lower cycle
         lowercyc = ceiling(unlistcycCT[[indij]]-(klag*mdim)-2)} #-2cycles }
       #list of lists of rss grey region
       rsslstarg[[i]][[j]] <- sum(lstarres[[i]][[j]]$residuals[lowercyc:uppercyc]^2) 
     }
   }
  ##rsslstarr THRESHOLD CT for +/- 2 cycles
   #add NA to those that are outside lower bound first
   unlistcycCTthr[which(unlistcycCTthr < (0.5+klag*mdim+2))] <- NA
   #rss for red box
   rsslstarr <- vector('list', sublength) 
   for(i in 1:(sublength)){
     for(j in 1:(replength/sublength)){
       #combined indicator with i and j
       indij = ((replength/sublength)*(i-1))+j  
       #klag-adjusted upper cycle b/c unlistcycCT[[indij]] is actual cycles for 3:46 (klag=2)
       #unlistcycCT[[indij]] - klag to fit lstarres$residuals 1:44 formmat
       #if(sum(is.na(unlistcycCTthr[[indij]])) > 0){
       #  rsslstarr[[i]][[j]] <- NA
       #  print(paste(k, "-", indij, "no CTth"))
       #}
       #filter out low values, change to NA unlistcycCTthr
       if(sum(isTRUE(unlistcycCTthr[[indij]] > cyclength[[i]]-(klag*mdim+1))) > 0){ #strict upper bound
         uppercyc = cyclength[[i]]-(klag*mdim)
         lowercyc = ceiling(unlistcycCTthr[[indij]]-(klag*mdim)-2)}
       else if(sum(is.na(unlistcycCTthr[[indij]])) == 0){ 
         uppercyc = floor(unlistcycCTthr[[indij]]-(klag*mdim)+2)  #+2 cycles
         #klag-adjusted lower cycle
         lowercyc = ceiling(unlistcycCTthr[[indij]]-(klag*mdim)-2)} #-2cycles 
       else{rsslstarr[[i]][[j]] <- NA}
       #list of lists of rss red region
       rsslstarr[[i]][[j]] <- sum(lstarres[[i]][[j]]$residuals[lowercyc:uppercyc]^2) 
     }
   }
   #unlisted 1:(replength) of rss grey region
   unlist.rsslstarg= unlist(rsslstarg, recursive = FALSE)
   unlist.rsslstarr= unlist(rsslstarr, recursive = FALSE)
   ##Creating RSS Matrix
   rss.mat <- matrix(c(unlist.rsslstar, unlist.rsslstarg, unlistcycCT, 
                                        unlist.rsslstarr, unlistcycCTthr), ncol=5)
  
 ###PART 2: LSTAR Parameter Coefficients  
   #generate coefficients of LSTAR model
   unlparams <- function(x){
     if(sum(is.na(x$fitted.values)) == 0){
       unlist(x$model.specific$par, recursive = FALSE)
     }
     else{ rep(NA, ((2*mdim)+4)) }
   }
   #list of replength=40 of coefficients
   lstarparamsl <- unlist(lapply(lstarres, function(x) lapply(x, unlparams)), recursive = FALSE)
   #matrix of lstar coefficients
   lstarparams <- do.call("rbind", lstarparamsl)
   ##LSTAR Model Parameters
   lstarparams.mat <- matrix(lstarparams, ncol=((2*mdim)+4)) 
   #klag*mdim for phiL and phiH, 4 from const.L, const.H, gamma, th
  
 ###PART 3: Durbin Watson Statistics
   #list of replength=40 of LSTAR models
   unlist.lstarmod <- unlist(lstarres, recursive=FALSE)
   cyclength <- lapply(subs, nrow)
   #replicated cycle length for each replicate
   repcyclength <- unlist(lapply(cyclength, function(x) rep(x, (replength/sublength))))
   #unlist.repcyc is the klag start of cycles
   unlist.repcyc <- list()
   for(i in 1:replength){
     unlist.repcyc[[i]] <- (1+klag*mdim):repcyclength[[i]] 
   }
   #output of durbin-watson for amplification and residuals
   reg.amp <- list() ; reg.res <- list() #dynlm amp and res
   for(i in 1:replength){
     if(sum(is.na(unlist.lstarmod[[i]]$fitted.values)) == 0){ #non-NAs (LSTAR fits)
       reg.amp[[i]] <- dynlm(as.numeric(unlist.lstarmod[[i]]$fitted.values) ~ unlist.repcyc[[i]]) #problems with zoo
       reg.res[[i]] <- dynlm(unlist.lstarmod[[i]]$residuals ~ unlist.repcyc[[i]])
     }
     else{
       reg.amp[[i]] <- NA
       reg.res[[i]] <- NA
     }
   }
   ##Finding Durbin Watson Statistics of Amplification and Residuals  
   durbinfunc <- function(x){
     tryCatch({
       durbinWatsonTest(x) #durbin watson for unlisted
     }, error=function(e) {
       return(list(r=NA, dw=NA, p=NA)) #error when durbinWatson(NA)
     })
   }
   ampdurbwat <- lapply(reg.amp, durbinfunc)
   resdurbwat <- lapply(reg.res, durbinfunc)
   #printing only autocorr, dw, and p-value
   matrixdw <- function(x){
     matrix(c(x$r, x$dw, x$p), ncol=3)
   }
   ##Creating DW-Statistics Matrix
   ampdw.mat <- do.call(rbind, lapply(ampdurbwat, matrixdw))
   resdw.mat <- do.call(rbind, lapply(resdurbwat, matrixdw))
   dw.mat <- cbind(ampdw.mat, resdw.mat)
  
 ###PART 4: Finding the Feature Set
   #FeatureSet for each group
   unlist.tst <- lapply(tst, function(x) unlist(x, recursive=FALSE)$FeatureSet) #rep in sub has same FeatureSet 
   # featsetl <- list()
   #  for(i in 1:10){
   #     featsetl[[i]] <- unlist.tst[[i]][seq(1, length(unlist.tst[[i]]), cyclength[[i]]*4)]
   #   }
   #final featureset vector
   #featset <- unlist(lapply(unlist.tst, function(x) rep(unique(x),(replength/sublength))))
  
  ##Rbind for Entire tst to make DF
   #rbindedtst <- list() 
   #for(i in 1:10){
   #   rbindedtst[[i]] <- do.call(rbind, tst[[i]])
   #  #rbindedtst[[i]] <- rbindedtst[[i]][order(rbindedtst[[i]]$SampleID),]
   #}
   #rbindedtst <- lapply(tst, function(x) do.call(rbind, x))
   #rbindedtst.df <- do.call(rbind, rbindedtst)
   #rbindedtst.df <- rbindedtst.df[order(rbindedtst.df$SampleID),]
   
   #cyclength = length of each subset ; cyclengtht = total length
   cyclengtht = lapply(cyclength, "*", (replength/sublength))
   #cumcyclength is cumulative across subsets
   cumcyclength = cumsum(cyclengtht) #accumulative 
   #unlisting tst
   tst.unl <- unlist(lapply(tst, function(x) unlist(x)))
   #selecting only FeatureSet
   featset <- tst.unl[grep("FeatureSet", names(tst.unl))]
   cycfeatset <- featset[cumcyclength]
   #each subset's featureset replicated
   featsetf <- list()
   for(i in 1:sublength){
     featsetf[[i]] <- rep(cycfeatset[[i]], (replength/sublength))
   }
  ##Finalized FeatureSet Matrix 
   featset.mat <- matrix(unlist(featsetf), ncol=1)
   
 ###PART 5: Finalizing Matrix Output   
  indk2 = replength*k ; indk1 = indk2-(replength-1)
  res[indk1:indk2, "FeatureSet"] <- featset.mat
  res[indk1:indk2, 5:10] <- dw.mat
  res[indk1:indk2, 11:13] <- rss.mat
  
  if(k==1){
  resp <- data.frame(lstarparams.mat)
  }
  else{
  resp[indk1:indk2, 1:(1+(2*mdim+4-1))] <- lstarparams.mat
  }
  colnames(resp) <- names(lstarparamsl[[1]])
  
    } #k files
  allres <- cbind(res, resp)
  return(allres)
}
