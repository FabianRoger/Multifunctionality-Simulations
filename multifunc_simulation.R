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

null.model <- data.frame(SP = as.factor(rep(LETTERS[1:10],each=10)),
                         Func = as.factor(rep(paste("F",formatC(c(1:10),width=2,flag="0"),sep="_"),10)),
                         Value = as.numeric(round(runif(100,0,1),digits=2)))
                         #Value = as.numeric(round(rnorm(100,2,1),digits=2)))

# function to calculate mixture values of function at given richness
avfunc <- function(R) {
  species.sample <- sample(levels(null.model$SP),R,replace = FALSE)
  sample.function.val <- null.model[null.model$SP %in%  species.sample,]
  average.function <- ddply(sample.function.val, .(Func), summarise, av_Func = mean(Value))
  return(average.function)
  }

# empty dataframe to store results from loop
mixture.null <- data.frame(Func = character(), 
                           av_Func = numeric(),
                           Richness = numeric(),
                           rep = numeric())

# populate dataframe for all richness levels from 1:10, with 100 replicates at each richness level
for( R in 1:10) {
  for( replicate in 1:100) {
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
  facet_wrap(~Func,nrow=2) +
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

#add on the new (standardized) functions along with the averaged multifunctional index
mixture_wide <- cbind(mixture_wide, getStdAndMeanFunctions(mixture_wide, colnames(mixture_wide)[3:12]))

#plot it
ggplot(mixture_wide, aes(x=Richness, y=meanFunction))+
    geom_point(size=3)+
    theme_bw(base_size=15)+
    stat_smooth(method="lm", colour="black", size=2) +
    xlab("\nSpecies Richness") +
    ylab("Average Value of Standardized Functions\n")

#reshape for plotting everything with ggplot2 
AvForPlot <- melt(mixture_wide[,c(1,13:23)], id.vars="Richness")

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

mixedThresh <- getFuncsMaxed(mixture_wide, colnames(mixture_wide)[3:12], threshmin=0.05, threshmax=0.99, 
                           prepend=c("Richness","replicate"), maxN=1)

gcPlot_mixed <- subset(mixedThresh, mixedThresh$thresholds %in% qw(0.2, 0.4, 0.6, 0.8))
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
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red")

# plot the slopes of the relationship against Threshold values

# note that I take fun = lm as glm doesn't converge for all thresholds over 67%
# however, the parts that converge look almost identical so it's no big deal

mixedLinearSlopes<-getCoefTab(funcMaxed ~ Richness, fun = lm,  data=mixedThresh, coefVar="Richness")

ggplot(mixedLinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate- 1.96*germanyLinearSlopes[["Std. Error"]],
                                 ymax=Estimate+1.96*germanyLinearSlopes[["Std. Error"]])) + 
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Change in Number of Functions per Addition of 1 Species\n") + xlab("\nThreshold (%)") +
  stat_abline(intercept=0, slope=0, lwd=1, linetype=2) + 
  theme_bw(base_size=14)

# highlight critical slopes

# We can't calculate Tmin, Tmax and Tmde because the calculation rely on the fact that the slopes get non-significant at
# some point. Becasue these are simulated data with 100 replicates for each richness level, this never happens (all slopes
# are significantly differnt from 0 at all threshold levels) . 




