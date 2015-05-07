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

#### define numper of Species (1 : 10000) ###

spec.num <- 5
### define number of Functions ##

func.num <- 3

########## define distribution of function values ########## 

## fixed sequence  ######## mean will depend on funcnum. have to fix that later
Funcval_seq <- c()
for (i in 1:spec.num) {
  Funcval_seq <-c(Funcval_seq, sample(seq(0.05,0.95,1/func.num),func.num,replace=F))
}

## random from uniform
Funcval_unif <- c()
for (i in 1:spec.num) {
  Funcval_unif <-c(Funcval_unif, runif(func.num,0,1))
}

## random from normal with mean of 4 and sd 1
Funcval_norm <- c()
for (i in 1:spec.num) {
  Funcval_norm <-c(Funcval_norm, rnorm(func.num,4,1))
}

## random from gamma 
Funcval_gamm <- c()
for (i in 1:spec.num) {
  Funcval_gamm <-c(Funcval_gamm, rgamma(func.num,1))
}

## random from beta 
Funcval_beta <- c()
for (i in 1:spec.num) {
  Funcval_beta <-c(Funcval_beta, rbeta(func.num,5,1))
}

## choose distribution

#Funcval <- Funcval_seq
Funcval <- Funcval_unif
#Funcval <- Funcval_norm
#Funcval <- Funcval_gamm
#Funcval <- Funcval_beta

######## German biodepth data #####

# load Biodepth data
data(all_biodepth)

# extract German data
Germany <- all_biodepth[all_biodepth$location == "Germany",]

# extract monocultures
Germany_mono <- Germany[Germany$Diversity == 1,]




## plot

hist(Funcval, breaks=10)

########## define null data set with parameters chosen above (Species, number functions, function values) #########

spec <- paste(LETTERS,rep(1:100,each=100),sep = "_")

null.model <- data.frame(SP = as.factor(rep(spec[1:spec.num],func.num)),
                         Func = as.factor(rep(paste("F",formatC(c(1:func.num),width = 2,flag = "0"),sep = "_"),each = spec.num)),
                         Value = as.numeric(round(Funcval,digits = 2)))
                         
# plot function correlation
null.wide <- dcast(null.model, SP ~ Func, value.var="Value")
ggpairs(null.wide[,-1], lower=list(continuous = "smooth"))

# calculate and plot possible combinations
possible.comb <- data.frame(Richness = rep(NA,spec.num), numb.comb =  rep(NA,spec.num))

for (p in 1: spec.num) {
  possible.comb$Richness[p] <- p
  possible.comb$numb.comb[p] <- prod(seq(spec.num, spec.num-(p-1),-1)) / prod(seq(1, spec.num-(spec.num-p),1))
}

ggplot(possible.comb, aes(x=Richness, y=numb.comb, label=numb.comb))+
  geom_point(size=3)+
  geom_text(hjust=1.5)+
  labs(y="number of possible species combinations", 
       title = paste("number of possible species combinations \n species pool = ", spec.num))+
  theme_bw(base_size=15)

######## function to calculate mixture values of function at given richness (unweighted average)
avfunc_unweighted <- function(R) {
  species.sample <- sample(levels(null.model$SP),R,replace = FALSE)
  sample.function.val <- null.model[null.model$SP %in%  species.sample,]
  average.function <- ddply(sample.function.val, .(Func), summarise, av_Func = mean(Value))
  return(average.function)
  }

######## function to calculate mixture values of function at given richness with selection effect

# I calculate the average function value in the mixture by weightin the species by the
# species' value of a specific function ( i.e. "biomass")

avfunc_selection <- function(R) {
  species.sample <- sample(levels(null.model$SP),R,replace = FALSE)
  sample.function.val <- null.model[null.model$SP %in%  species.sample,]
  Weights <- sample.function.val[sample.function.val$Func == "F_01",]$Value
  sample.function.val$weights <- rep( c( Weights / sum( Weights)), func.num)
  average.function <- ddply(sample.function.val, .(Func), summarise, av_Func = weighted.mean(Value,weights))
  return(average.function)
}

######## function to calculate mixture values of function at given richness with complementarity effect
avfunc_complementarity <- function(R) {
  species.sample <- sample(levels(null.model$SP),R,replace = FALSE)
  sample.function.val <- null.model[null.model$SP %in%  species.sample,]
  sample.function.val$Value <- sample.function.val$Value * 1.01^R
  average.function <- ddply(sample.function.val, .(Func), summarise, av_Func = mean(Value))
  return(average.function)
}

########## choose function #########

avfunc <- avfunc_unweighted
#avfunc <- avfunc_selection
#avfunc <- avfunc_complementarity


# empty dataframe to store results from loop
mixture.null <- data.frame(Func = character(), 
                           av_Func = numeric(),
                           Richness = numeric(),
                           rep = numeric())



# create species matrix with all Diversity levels, enevn replication and 

# all monocultures
# total number of plots 
# each monoculture, at each richness level, each species is present in at least 3 plots in random combinations
# each richness level hass all possible combination if < 50, min 50 and half of all possible combinations if < 100
# and maximum 200 combinations



nplot <- spec.numb

Spec.mat <- matrix( data = 0, nrow = nplot, ncol = spec.num, dimnames = list( c(1:nplot), levels(null.model$SP)))

possible.comb$nplot  <- possible.comb$numb.comb * 

  spec.num
















# populate dataframe for all richness levels from 1:10, with 100 replicates at each richness level
for( R in 1:spec.num) {
  for( replicate in 1:10) {
    temp.func <- avfunc(R)
    temp.func$Richness <- R
    temp.func$replicate <- replicate
    mixture.null <- rbind(mixture.null,temp.func)
  }
}

# cast dataframe in wide format to comply with the dataformat used in the multifunc package
mixture_wide <- dcast( mixture.null, Richness + replicate ~ Func, value.var = "av_Func")

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

# i dont calculate slopes and p-values as I'm primarily intrested in the graph here

ggplot(mixture.null, aes(x = Richness, y = av_Func))+
  geom_point(size=3)+
  facet_wrap(~Func) +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) +
  xlab("\nSpecies Richness") +
  ylab("Value of Function\n") +
  theme(panel.grid = element_blank())
  

#####################
# OVERLAPP APPROACH #
#####################

#(I jump over this one for now)

######################
# AVERAGING APPROACH #
######################

# extract function names
func.names <- levels(mixture.null$Func)

# add on the new (standardized) functions along with the averaged multifunctional index
mixture_wide <- cbind(mixture_wide, getStdAndMeanFunctions(mixture_wide, func.names))

#plot it
ggplot(mixture_wide, aes(x=Richness, y=meanFunction))+
    geom_point(size=3)+
    theme_bw(base_size=15)+
    stat_smooth(method="lm", colour="black", size=2) +
    xlab("\nSpecies Richness") +
    ylab("Average Value of Standardized Functions\n")

#reshape for plotting everything with ggplot2 
AvForPlot <- melt(mixture_wide[,c("Richness",paste(func.names,"std",sep="."))], id.vars="Richness")

#plot it 
ggplot(AvForPlot, aes(x=Richness,y=value))+
  geom_point(size=3)+ 
  facet_wrap(~variable,nrow=2) +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) + 
  xlab("\nSpecies Richness") +
  ylab("Standardized Value of Function\n")

######################
# THRESHOLD APPROACH #
######################

mixedThresh <- getFuncsMaxed(mixture_wide, func.names, threshmin=0.05, threshmax=0.99, 
                           prepend=c("Richness","replicate"), maxN=1)

gcPlot_mixed <- subset(mixedThresh, mixedThresh$thresholds %in% qw(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
gcPlot_mixed$percent <- paste(100*gcPlot_mixed$thresholds, "%", sep="")

# note: Jerret uses GLMs to specify the best fit lines. As they don't converge in all cases for this dataframe,
# I use just lm(). For the cases where the glm converged, the differences seem to be minimal.

# plot 4 selected threshold values
ggplot(gcPlot_mixed, aes(x=Richness, y=funcMaxed))+
  geom_point()+
  facet_wrap(~percent)+
  stat_smooth(method="lm",colour="red", lwd=1.2 )+
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

mixedLinearSlopes<-getCoefTab(funcMaxed ~ Richness, fun = lm,  data=mixedThresh, coefVar="Richness")

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



