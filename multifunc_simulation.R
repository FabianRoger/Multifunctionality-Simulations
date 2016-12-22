#library(devtools)
#install_github("multifunc", "jebyrnes")

library(plyr)
library(multifunc)
library(reshape2)
library(ggplot2)
library(GGally)
library(gridExtra)

##############################################################################

# simulation for multifunc project. 2015-05-04 #

##############################################################################

# create dataset with 10 species and 10 functions with mixtures

# each species has a random function value for each of the function
# scaled between 0 and 1; values drawn form uniform distribution #dunif

# null model: substitutive and purely additive. multifunnc is average of weighted
# average of single func (here: completely even so just arithmetic mean)

########### load functions #############

source("functions.R")

####################


#### define numper of Species (1 : 12) ###

specnum <- 11
### define number of Functions ##

funcnum <- 5

# set up species matrix with specnum species 

SpecMat <- SpeciesMatrix(specnum = specnum)


# set up functio matrix with specnum Species, funcnum functions
# the functino values are drawn from the distribution specified in distribution
# (all distribution listed in ?Distributions should work), all paramters
# besides the number of samples to draw from the distribution must be specified

FuncMat <- FunctionValue(specnum = specnum, funcnum = funcnum, spec = NULL, 
                         func = NULL, distribution = "runif", 0, 1)

# calculate plotwise function averages
# three methods can be specified: average, complementarity and selection
# complementarity and selection have further details, see the documentation in 
# the function code for details. 

# example with selection effect
AvFunc <- AverageFunction(SpecMat, FuncMat, method = "av", selfac = 2 )

# example with complementarity effect
AvFunc <- AverageFunction(SpecMat, FuncMat, method = "com", comp = 1.05, compfunc = "Func_01" )


######## Biodepth data  #####

# load Biodepth data
data(all_biodepth)

# function consodered in JByrnes et al 
allVars<-qw(biomassY3, root3, N.g.m2, N.Soil, cotton3) 
varIdx<-which(names(all_biodepth) %in% allVars)

# format Biodepth data to match AvFunc structure
Biodepth <- cbind( all_biodepth[ , 26:ncol( all_biodepth)],
                   all_biodepth[ , c(1,5,8)],
                   all_biodepth[ , varIdx])

colnames(Biodepth)[which(colnames(Biodepth) == "Diversity")]  <- "Richness"

###### choose country ########
AvFunc <- Biodepth[Biodepth$location == "Germany",][,which(
  colnames(Biodepth) != "location")]

specnum <- ncol(all_biodepth) - 25


#################################
# null model with Biodepth data #

#################################

levels(all_biodepth$location)

# select country
biodepth.null <- all_biodepth[ all_biodepth$location == "Germany" , ]

# extract monocultures
monocultures <- biodepth.null[biodepth.null$Diversity == 1, ]

# which species present in this country 
spec <- names(which(colSums(monocultures[,26:ncol(monocultures)]) != 0))
specnum <- length(spec)

# select function
func <- c("biomassY3","root3","N.g.m2","N.Soil","cotton3" )

funcNA <- which(is.na(monocultures[,func]), arr.ind=T)

if (nrow(funcNA) > 0) {
func <- func[ -c(unique(funcNA[,2]))]
}

funcnum <- length(func)

# subset moncultures for reelvant columns
monocultures <- monocultures[,c(spec,func)]

# extract monoculture values for all species and shape to long dataframe
 # comnpatible with AverageFunction()

monocultures$Species <- NA

for ( i in 1:nrow(monocultures)) {
  W <- which(monocultures[i,spec] == 1)
  monocultures$Species[i]  <- colnames(monocultures)[W]
}

function.val <- melt(monocultures[,! names(monocultures) %in% spec], 
                  id.var = "Species")

colnames(function.val) <- c("Species", "Functions", "Funcval")

# standardize function values between 0 and 1
function.val$Funcval <- unlist(by(function.val, function.val$Functions, 
                                  function(x) {x$Funcval <- x$Funcval/max(
                                    x$Funcval)}))


# create species matrix 
Species.Matrix <- SpeciesMatrix( spec = spec)

# add average fucntion
AvFunc <- AverageFunction( Species.Matrix, function.val, method = "av",
                           comp = 1.05, selfac = 1, selfunc = "biomassY3", 
                           compfunc = "all")


#################
# In the following I replicate the analysis steps described 
# In:
# Supplementary Information 1: Using the multifunc package for analysis of 
# Biodiversity Ecosystem Multifunctionality Relationships in R
# from: 
# Byrnes et al 2014; Investigating the relationship between biodiversity and 
# ecosystem multifunctionality: challenges and solutions

# The analysis above uses the german BIODEPTH data. 
# For comarision, I replicate all performed here also with the original data
# see Biodepth_example_JByrnes.R

# If I have to deviate from the analysis I specify it in the script
###################

############################
# SINGLE FUNCTION APPROACH #
############################

# reshape for plotting

AvFunc_long <- melt(AvFunc[,-c(1:specnum)], id.vars=c("plot", "Richness"))

ggplot(AvFunc_long, aes(x = Richness, y = value))+
  geom_point(size=3)+
  facet_wrap(~variable) +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) +
  xlab("\nSpecies Richness") +
  ylab("Value of Function\n") +
  theme(panel.grid = element_blank())
  

######################
# AVERAGING APPROACH #
######################

# extract function names
func.names <- colnames( AvFunc[ ( specnum + 3) : ncol( AvFunc)])

# add on the new (standardized) functions along with the averaged multifunctional index
AvFunc <- cbind(AvFunc, getStdAndMeanFunctions(AvFunc, func.names))

#plot it
ggplot(AvFunc, aes(x=Richness, y=meanFunction))+
    geom_point(size=3)+
    theme_bw(base_size=15)+
    stat_smooth(method="lm", colour="black", size=2) +
    xlab("\nSpecies Richness") +
    ylab("Average Value of Standardized Functions\n")

#reshape for plotting everything with ggplot2 
AvForPlot <- melt(AvFunc[,-c((specnum+3) : (specnum+funcnum+2))],
                  id.vars=colnames(AvFunc[1:(specnum+2)]))


# plot each function
ggplot(AvForPlot, aes(x=Richness, y=value))+
  geom_point(size=3)+ 
  facet_wrap(~variable,nrow=2) +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) + 
  labs( x = "\nSpecies Richness", y = "Standardized Value of Function\n")

#plot it (with colour code for each species) (comment out qw funciton to execute)

qw(
  
for (SP in spec) {
  gg <- ggplot(AvForPlot, aes_string(x="Richness", y="value", 
                                     colour = SP))+
    geom_point(size=3)+ 
    facet_wrap(~variable,nrow=2) +
    theme_bw(base_size=15)+
    stat_smooth(method="lm", colour="black", size=2) + 
    labs( x = "\nSpecies Richness", y = "Standardized Value of Function\n",
          title = SP)
  
  print(gg) 
  }

)


######################
# THRESHOLD APPROACH #
######################

mixedThresh <- getFuncsMaxed(AvFunc, func.names, threshmin=0.05, threshmax=0.99, 
                           prepend=c("plot","Richness"), maxN=7)

gcPlot_mixed <- subset(mixedThresh, mixedThresh$thresholds %in% qw(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
gcPlot_mixed$percent <- paste(100*gcPlot_mixed$thresholds, "%", sep="")

# note: Jerret uses GLMs to specify the best fit lines. As they don't converge in all cases for this dataframe,
# I use just lm(). For the cases where the glm converged, the differences seem to be minimal.

# plot 4 selected threshold values
ggplot(gcPlot_mixed, aes(x=Richness, y=funcMaxed))+
  geom_point()+
  facet_wrap(~percent)+
  stat_smooth(method="glm",colour="red", lwd=1.2 )+
  labs(y = expression("Number of Functions" >= Threshold), x = ("Species Richness"))+
  theme_bw(base_size=15)


# plot N Function over Threshold ~ Richness Slopes for all Thresholds
mixedThresh$percent <- 100*mixedThresh$thresholds 

ggplot(data=mixedThresh, aes(x=Richness, y=funcMaxed, group=percent)) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("Species Richness") +
  stat_smooth(method="lm", lwd=0.8, fill=NA, aes(color=percent))+
  theme_bw(base_size=14) +
  scale_y_continuous(limits=c(0,max(mixedThresh$funcMaxed)))+
  geom_hline(aes(yintercept = max(mixedThresh$funcMaxed)),size=1)+
  geom_hline(aes(yintercept = 0), size=1)+
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red")

# plot the slopes of the relationship against Threshold values

# note that I take fun = lm as glm doesn't converge for all thresholds over 67%
# however, the parts that converge look almost identical so it's no big deal

mixedLinearSlopes<-getCoefTab(funcMaxed ~ Richness, fun = lm,  data=mixedThresh, 
                               coefVar="Richness")

colnames(mixedLinearSlopes) <- c("thresholds", "Estimate",  "Std. Error", "t value", "Pr(>|t|)")


ggplot(mixedLinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate- 1.96*mixedLinearSlopes[["Std. Error"]],
                                 ymax=Estimate+1.96*mixedLinearSlopes[["Std. Error"]])) + 
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Change in Number of Functions per Addition of 1 Species\n") + xlab("\nThreshold (%)") +
  stat_abline(intercept=0, slope=0, lwd=1, linetype=2) + 
  theme_bw(base_size=14)

# highlight critical slopes

# We can't calculate Tmin, Tmax and Tmde because the calculation rely on the fact that the slopes get non-significant at
# some point. Becasue these are simulated data with 100 replicates for each richness level, this never happens (all slopes
# are significantly differnt from 0 at all threshold levels) . 



