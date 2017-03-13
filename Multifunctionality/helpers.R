######## function definitions #############


PossibleCombinations  <- function(specnum, maxrep) {
  # Calculates all possible combinations at each richness level 
  # from 1 to spec.num
  # Calculates suggested replication at each richness level as follows:
  # if possible combinations <= maxrep; suggested replication = all combinations
  # if possible combinations > maxrep; suggested maxrep replications
  # Calculates cumulative plotnumber
  #
  # Args:
  #   spec.num: integer, positive number giving the number of species in the
  #             species pool
  #   maxrep:   integer, positive number giving the maximum numer of replicated communities per richness level
  #
  # Returns:
  #   data.frame with four columns
  #   $Richness: richness levels from 1:specnum
  #   $numb.comb: number of possible combination at each Richness level
  #   $nrep: suggested replication at each richness level
  #   $sumrep: cummulative number of plots
  
  possible.comb <- data.frame(Richness = 1:specnum)
  possible.comb$numb.com <- choose(specnum, possible.comb$Richness)
  possible.comb$nrep <- possible.comb$numb.com
  if (max( possible.comb$nrep) > maxrep) {
    possible.comb[possible.comb$nrep > maxrep,]$nrep  <- maxrep
  }
  possible.comb$sumrep <- cumsum(possible.comb$nrep)
  
  return(possible.comb)
  
}

########################################################################################

SpeciesList <- function(specnum) {
  # creates character vector with specnum species names. If specnum ≤ 26 Species are 
  # named as A - Z, if specnum > 26, species are named as follows: 
  # A_01, B_01 (...) Z_01, A_02, B_02 (...)
  #
  # Args:
  #   specnum: integer, positive number giving the number of species to be named
  #
  # Returns: character vector with specnum species names
  
  L <- length(LETTERS)
  NUM <- ceiling(specnum/L)
  if (specnum <= L) {
    spec.list <- LETTERS[1:specnum]
  } else {
    spec.list <- paste(LETTERS, rep(formatC(1:NUM, width=2, flag="0"), each=L),
                       sep = "_")
    spec.list <- spec.list[1:specnum]
  }
  
  return(spec.list)
  
}

########################################################################################

FunctionList <- function(funcnum) {
  # creates character vector with funcnam function names. Function are named as 
  # follows: F 1, F 2 (...) 
  #
  # Args:
  #   funcnum: integer, positive number giving the number of Functions to be named
  #
  # Returns: character vector with funcnam Function names
  
  func.list <- paste("F", c(1:funcnum),sep = " ")
  
  return(func.list)
  
}

########################################################################################

SpeciesMatrix <- function(specnum, maxrep = 50) {
  # creates a species matrix with specnum species.
  # At each richness level, the replicates communities are all unique
  # The species names are generated with SpeciesList() 
  # and the experimental design is generated with the function PossibleCombinations()
  #
  #
  # Args:
  #   specnum: integer, positive number giving the number of species in the
  #             species pool
  #   maxrep: integer, positive number giving the maximum numer of replicated communities per richness level
  #
  # Returns:
  #   Species matrix with with presence absence data for species in plots
  #   rows are plots, species are column. Columnnames are species names, taken
  #   from specnum. Number of richness levels and number of replication
  #   per richness levels are taken from expdesign 
  
  spec <- SpeciesList(specnum) #list of species names
  
  expdesign <- PossibleCombinations( specnum, maxrep) #experimental design (number of replication per richness level)
  
  nplot <- sum(expdesign$nrep) #total number of replicates
  
  # empty species matrix
  Spec.mat <- matrix( data = 0, nrow = nplot, ncol = specnum, 
                      dimnames = list( c(1:nplot),spec))
  
  for (i in seq_along(expdesign$Richness)) { #fill species matrix with unique species combinations
    
    if (i == 1) { # fill monocultures
      for (n in 1:specnum) {
        Spec.mat[n,n]  <- 1
        
      }} else if (i < max(expdesign$Richness)) { # fill all other richness levels but highest
        
        
        M <- matrix( data = sample(spec, i), nrow = i, ncol = expdesign$nrep[i])
        M <- apply(M, 2, sort)
        
        for (n in 2: ncol(M)) { # create matrix with 
          
          next.spec <- sort(sample(spec, i))
          
          while( i %in% colSums(as.matrix(next.spec == M[, 1:(n-1)]))) {
            
            next.spec <- sort(sample(spec, i))
            
          }
          
          M[,n]<- next.spec
          
        }
        
        row.seq <- (expdesign$sumrep[i-1]+1) : expdesign$sumrep[i]
        
       
        for(k in row.seq) {
          Spec.mat[k, which(spec %in% M[, (k+1) - min(row.seq)])] <- 1
        }
        
      } else {Spec.mat[nrow(Spec.mat), ] <- 1}
  }
        
  
  return(Spec.mat)
  
  }

##################################################

FunctionValue <- function(specnum = NULL, funcnum = NULL, 
                          distribution = NULL, ...) {
  # function that assigns N species function values for M functions
  #
  # Args:
  #   specnum: number of species (numeric, integer)
  #   funcnum: number of function
  #   distribution: any of the distribution functions in the stats package
  #                 given as character string (e.g. "rnorm"). 
  #   ... further arguments to be passed to the distribution function
  
  # error handling
  if (! distribution %in% ls("package:stats")) {
    stop(" distribution must be one of the functions in the stats package.
         see ?Distributions for options")
  }
  
  distf <- get(distribution)
  
  Try <- try(distf(5,...))
  
  if(class(Try) == "try-error") {
    stop ( paste("distribution parameters seem to be wrongly speciefied see ?",
                 distribution," for help", sep=""))}
  ###
  
    spec <- SpeciesList( specnum)

    func <- FunctionList( funcnum)
  
  SpecFunc <- data.frame(Species = rep(spec,funcnum), 
                         Functions = rep(func, each = specnum), 
                         Funcval = NA)
  
  for (i in func) {
    SpecFunc[SpecFunc$Functions == i,]$Funcval <- distf(specnum,...)
  }
  
  if(TRUE %in% c(SpecFunc$Funcval < 0)) {
    warning("negative function values were generated wherefore all values were shifted by abs(min(function value))")
    SpecFunc$Funcval <- SpecFunc$Funcval + abs(min(SpecFunc$Funcval))
  }
  
  return(SpecFunc)
  
  }


##################################################

AverageFunction <- function(SPM, FUNC, method = "average", selfunc = "Func_01", selfac = 1,
                            compfunc = "all", CF = 2, r = 0.25) {
  # function to calculate the avreage function value of species
  # mixtures. 
  #
  # Args:
  #   SPM: a species matrix as prodcued by SpeciesMatrix(); if matrix is
  #        not presence/absence, it is transformed to p/a internally
  #   FUNC: a dataframe as produced by FunctionValue() giving the function 
  #         values for each species in long format. Species names must match 
  #         colnames of SPM. 
  #   method: one of "average", "complementarity" or "selection", specifying 
  #           how the mixture function should be calculated
  #   CF: "complementarity factor" for method "complementarity". each function
  #        value will be multiplied by CF * ( 1 - ( 1 - 1/CF ) * exp(1-Richness^r))
  #        the function growth to CF as Richness -> inf
  #   r: growthrate of amplification factor (how fast it growth to CF)
  #   selfunc: "selection function", function whichs function values are taken 
  #             to calculate the weighted mean for all functions
  #   selfac: "selection factor": factor by which the weights will be amplified 
  #   compfunc: if complementarity should only happen for some functions, the 
  #             functions can be specified here. default is all. 
  #
  # Returns:
  #   dataframe with Richness of each plot and average function value for each 
  #   function
  
  # error handling
  if (class(SPM) != "matrix" | is.numeric(SPM) != TRUE) {
    stop("SPM must be a numeric matrix")
  }
  
  if (NA %in% match(colnames(SPM), unique(FUNC$Species))) {
    stop(" not all species in species matrix represented in function dataframe")
  }
  
  spec <- colnames(SPM)
  
  specnum <- length(spec)
  
  # define function to be applied to each plot for method average or complementarity
  FUNC <- FUNC %>% spread(Functions, Funcval)
  
  if (is.na(pmatch(method, c("average", "complementarity"))) == FALSE) {
    MeanFunc <- function(x) {
      colMeans(FUNC[FUNC$Species %in% spec[which(x == 1)], ][,-1])
    } 
  }
  
  # define function to be applied to each plot for selection effect
  if (is.na(pmatch(method, "selection")) == FALSE) {
    
    spec.weights <- FUNC[, which(colnames(FUNC) == selfunc)]^selfac
    
    MeanFunc <- function(x) {
      Weights <- spec.weights[which(x == 1)]
      Weights <- Weights/sum(Weights)
      colSums(FUNC[FUNC$Species %in% spec[which(x == 1)], ][,-1]*Weights)
    } 
  }
  
  
  mean.functions <- t(apply(SPM, 1, MeanFunc))
  
  SPM <- as.data.frame(SPM)
  SPM$Richness <- rowSums(SPM)
  SPM <- cbind(SPM , mean.functions)
  
  # multiply with complementraity factor if method = "complementarity"
  
  if (is.na(pmatch(method, "complementarity")) == FALSE) {
    
    if ("all" %in% compfunc) {
      
      SPM[ , ( specnum + 2) : ncol( SPM)] <- SPM[ , ( specnum + 2) : ncol( SPM)] * 
        (CF * ( 1 - ( 1 - 1/CF ) * exp(1-SPM$Richness^r)))
      
    } else { SPM[ , which( colnames( SPM) %in% compfunc)] <- SPM[ , which( colnames( SPM) %in% compfunc)] *
      (CF * ( 1 - ( 1 - 1/CF ) * exp(1-SPM$Richness^r)))}
  }
  
  return(SPM)
  
}


##################################################



SlopeSummary <- function(MT) {
  # function to find at which threshold the estimated slope
  # reaches it's minimum
  # reaches it's maximum 
  # changes from positive to negative
  #
  # Args:
  #   MT: a dataframe as produced by multifunc::getCoefTab that needs to contain at least three columns, 
  #       one called "Estimate" that is taken as the column containing the slope estimates for each threshold,
  #       one called "Std. Error" that is taken to contain the corresponding standard errors and one called
  #       "thresholds" that contains the threshold levels. 
  #   
  #
  # Returns:
  #   named list containing the threshold level for which the sign change occours for both the estimate and the minimum estimate
  #   the threshold is always the threshold before which the sign change has occured
  
  
  ### error handling
  if (class(MT) != "data.frame") {
    stop("the function expects a dataframe containing at least the following columns: 
         `Estimate`, `Std. Error` and `thresholds`")
  }
  
  
  if (NA %in% match(c("Estimate", "Std. Error",  "thresholds"), colnames(MT))) {
    stop("one or more of the columns `Estimate`, `Std. Error` and `thresholds` are missing form the dataframe")
  }
  ###
  
  #MT$min.est <- MT$Estimate - 1.96*MT$`Std. Error`
  MT$est_sign_change <- "no"
  #MT$min.est_sign_change <- "no"
  
  for (n in 2: nrow(MT) ){
    ifelse( sum(MT$Estimate[ c( n-1, n)] >= 0) == 1, 
            MT$est_sign_change[n] <- "yes",
            MT$est_sign_change[n] <- "no")
    
    #ifelse( sum(MT$min.est[ c( n-1, n)] >= 0) == 1, 
     #       MT$min.est_sign_change[n] <- "yes",
      #      MT$min.est_sign_change[n] <- "no")
  }
  
  MAX <- MT[ which(MT$Estimate == max(MT$Estimate)), ]$thresholds
  MIN <- MT[ which(MT$Estimate == min(MT$Estimate)), ]$thresholds
  
  
  Signchange.R <- list(maximum = MAX,
                       minimum = MIN,
                       Estimate_sign_change = MT[ MT$est_sign_change == "yes", ]$thresholds)#,
                     #min.Estimate_sign_change = MT[ MT$min.est_sign_change == "yes", ]$thresholds)
  
  if(length(Signchange.R$Estimate_sign_change) == 0) {Signchange.R$Estimate_sign_change <- NA}
 # if(length(Signchange.R$min.Estimate_sign_change) == 0) {Signchange.R$min.Estimate_sign_change <- NA}
  
  return(Signchange.R)
}


########### function taken form the multifunc package ###################

# The functions below are taken from the multifunc package and written by Jarret Byrnes (https://github.com/jebyrnes/multifunc)

# They are included here in a modified form as 
# 1) the package is not available on CRAN for the moment
# 2) the package loads plyr which conflicts with dplyr if loaded afterwards

#full citation:

# Byrnes, J. E. K., Gamfeldt, L., Isbell, F., Lefcheck, J. S.,
# Griffin, J. N., Hector, A., Cardinale, B. J., Hooper, D. U., Dee, L. E.,
# Emmett Duffy, J. (2014), 
# Investigating the relationship between biodiversity and ecosystem multifunctionality: 
# challenges and solutions. 
# Methods in Ecology and Evolution, 5: 111–124. doi: 10.1111/2041-210X.12143

# note that I re-wrote the functions below in order to not rely on plyr

#getFuncsMaxed
############################################

getFuncsMaxed<-function(adf, vars=NA, threshmin=0.05, threshmax=0.99, threshstep=0.01, proportion=F, prepend="Diversity", maxN=1) {
  
  thresholds <- seq(threshmin,threshmax,threshstep)
  
  ret_list <- lapply(thresholds, function(x) {
    
    df <- getFuncMaxed(adf,
                       vars=vars,
                       thresh= x, #the threshold
                       #proportion=proportion,
                       maxN=maxN,
                       prepend=c(prepend))
    
    df$thresh <- as.numeric(x)
    
    #to match original output exactly
    df$thresholds <- df$thresh
    df <- df[,c(5,1:4)]
    
    return(df)
  })
  
  ret <- bind_rows(ret_list)
  
}



#getFuncMaxed
################
getFuncMaxed<-function(adf, vars=NA, thresh=0.7, prepend="Diversity", maxN=1){
  if(is.na(vars)[1]) stop("You need to specify some response variable names")
  
  #scan across all functions, see which are >= a threshold
  #funcMaxed<-rowSums(colwise(function(x) x >= (thresh*max(x, na.rm=T)))(adf[,which(names(adf)%in%vars)]))
  
  getMaxValue<-function(x){
    l<-length(x)    
    mean( sort(x, na.last=F)[l:(l-maxN+1)], na.rm=T)
  }
  
  funcMaxed<-rowSums(apply(adf[,which(names(adf)%in%vars)], 2, function(x) x >= thresh*getMaxValue(x)))
  
  #bind together the prepend columns and the functions at or above a threshold
  ret<-data.frame(cbind(adf[,which(names(adf) %in% prepend)], funcMaxed))
  names(ret) <- c(names(adf)[which(names(adf) %in% prepend)], "funcMaxed")
  
  #how many functions were considered
  ret$nFunc<-length(vars)
  
  ret
}
############################################



# getCoefTab
############################################
getCoefTab<-function(eqn, fun=glm, data, groupVar="thresholds", coefVar, ...){
  
  getCoef <- function(adf, ...) {
    
    if(length(unique(adf[[ as.character(eqn[[2]]) ]]))==1){ #in case all functions perform exactly the same, the slope is a flat line
      
      df <- data.frame(t(c(0,0,NA,1)))
      colnames(df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      return(df)
      
    } 
    
    options(warn=2)
    
    #use try-catch in case there are errors
    aFit<-try(fun(eqn, data=adf, ...))
    
    options(warn=0)
    
    #if there was a problem, catch it and just return NAs for this coefficient
    if("try-error" %in% class(aFit)) return(rep(NA, 4))
    
    if("glm" %in% class(aFit)) {    if(!aFit$converged) return(rep(NA, 4))	}
    
    coefInfo<-summary(aFit)$coef
    
    idx<-which(rownames(coefInfo) == coefVar)
    
    df <- data.frame(t(coefInfo[idx,]))
    colnames(df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    
    return(df)
  } 
  
  thresh_list <- split(data, data[, groupVar])
  
  ret_list <- lapply(thresh_list, getCoef)
  ret <- do.call("rbind", ret_list)
  
  ret$thresholds <- as.numeric(names(thresh_list))
  ret <- ret[, c("thresholds", "Estimate",  "Std. Error", "t value", "Pr(>|t|)")]
  
  return(ret)
  
}
############################################

  

