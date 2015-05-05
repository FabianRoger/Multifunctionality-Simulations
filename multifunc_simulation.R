#library(devtools)
#install_github("multifunc", "jebyrnes")

library(plyr)
library(multifunc)
library(reshape2)
library(ggplot2)
library(gridExtra)

##############################################################################

# simulation for multifunc project. 2015-05-04 #

##############################################################################

# create dataset with 10 species and 10 functions with mixtures

# each species has a random function value for each of the function
# scaled between 0 and 1; values drawn form uniform distribution #dunif

# null model: substitutive and purely additive. multifunnc is average of weighted
# average of single func (here: completely even so just arithmetic mean)

#### define numper of Species (1 : 10000) ###

spec.num <- 20

### define number of Functions ##

func.num <- 15

########## define distribution of function values ########## 

## fixed sequence  ######## mean will depend on funcnum. have to fix that later
Funcval_seq <- c()
for (i in 1:10) {
  Funcval_seq <-c(Funcval_seq, sample(seq(0.05,0.95,1/func.num),10,replace=F))
}

## random from uniform
Funcval_unif <- c()
for (i in 1:10) {
  Funcval_unif <-c(Funcval_unif, runif(func.num,0,1))
}

## random from normal with mean of 4 and sd 1
Funcval_norm <- c()
for (i in 1:10) {
  Funcval_norm <-c(Funcval_norm, rnorm(func.num,4,1))
}

## random from gamma 
Funcval_gamm <- c()
for (i in 1:10) {
  Funcval_gamm <-c(Funcval_gamm, rgamma(func.num,1))
}

## choose distribution

#Funcval <- Funcval_seq
#Funcval <- Funcval_unif
#Funcval <- Funcval_norm
Funcval <- Funcval_gamm


########## define null data set with parameters chosen above (Species, number functions, function values) #########

  