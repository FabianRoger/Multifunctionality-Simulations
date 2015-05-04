library(multifunc)
library(ggplot2) 
library(gridExtra)

################################################################################
# analysis from JByrnes et al for comparision with modeled multifunc data

## Supplementary Information 1: Using the multifunc package for analysis of 
# Biodiversity Ecosystem Multifunctionality Relationships in R
# from: 
# Byrnes et al 2014; Investigating the relationship between biodiversity and 
# ecosystem multifunctionality: challenges and solutions

# The analysis below uses the german BIODEPTH data. 
################################################################################

# load and subset data

data(all_biodepth)
allVars<-qw(biomassY3, root3, N.g.m2, light3, N.Soil, wood3, cotton3) 
varIdx<-which(names(all_biodepth) %in% allVars)


#Now, specify what data we're working with - Germany
#and the relevant variables we'll be working with

germany<-subset(all_biodepth, all_biodepth$location=="Germany")
vars<-whichVars(germany, allVars) 
species<-relevantSp(germany,26:ncol(germany))
spIDX <- which(names(germany) %in% species) #in case we need these

############################
# SINGLE FUNCTION APPROACH #
############################

germanyForPlotting<-melt(germany[,c(8,which(names(germany) %in% vars))], id.vars="Diversity")
germanyForPlotting$variable <- factor(germanyForPlotting$variable)

#make the levels of the functions into something nice for plots 
levels(germanyForPlotting$variable) <- c('Aboveground Biomass', 'Root Biomass', 'Cotton Decomposition', 'Soil Nitrogen', 'Plant Nitrogen')


germanyFits <- dlply(germanyForPlotting, .(variable), function(x) lm(value ~ Diversity, data=x))

germanyLabels <- data.frame(variable = levels(germanyForPlotting$variable), 
                            Diversity=7,
                            value=c(2000, 1200, 1.5, 6, 20),
                            lab=paste(letters[1:5], ")", sep=""),
                            r2 = sapply(germanyFits, function(x) summary(x)$r.squared), 
                            p = sapply(germanyFits, function(x) anova(x)[1,5]))
                           
germanyLabels$labels <- with(germanyLabels, paste(lab, "p =", round(p,3),expression(R^2), "=", round(r2,2), sep=" "))
germanyLabels$labels <- gsub("p = 0 ", "p < 0.001 ", germanyLabels$labels)

ggplot(germanyForPlotting, aes(x=Diversity, y=value)) + 
  geom_point(size=3)+
  facet_wrap(~variable, scales="free") +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) +
  xlab("\nSpecies Richness") +
  ylab("Value of Function\n") +
  geom_text(data=germanyLabels, aes(label=labels)) + 
  theme(panel.grid = element_blank())

#####################
# OVERLAPP APPROACH #
#####################

#(I jump over this one for now)


######################
# AVERAGING APPROACH #
######################

#add on the new (standardized) functions along with the averaged multifunctional index
germany<-cbind(germany, getStdAndMeanFunctions(germany, vars)) 


#plot it
ggplot(germany, aes(x=Diversity, y=meanFunction))+
  geom_point(size=3)+
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) +
  xlab("\nSpecies Richness") +
  ylab("Average Value of Standardized Functions\n")

#reshape for plotting everything with ggplot2 
germanyMeanForPlotting<-melt(germany[,c(8,129:134)], id.vars="Diversity")

#nice names for plotting
levels(germanyMeanForPlotting$variable) <- c('Aboveground Biomass', 'Root Biomass', 'Cotton Decomposition', 'Soil Nitrogen', 
                                             'Plant Nitrogen', 'Mean Multifuncion Index')
#plot it 
ggplot(aes(x=Diversity,y=value),data=germanyMeanForPlotting)+
  geom_point(size=3)+ 
  facet_grid(~variable) +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) + xlab("\nSpecies Richness") +
  ylab("Standardized Value of Function\n")


######################
# THRESHOLD APPROACH #
######################

germanyThresh<-getFuncsMaxed(germany, vars, threshmin=0.05, threshmax=0.99, prepend=c("plot","Diversity"), maxN=7)
gcPlot<-subset(germanyThresh, germanyThresh$thresholds %in% qw(0.2, 0.4, 0.6, 0.8))
gcPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")


# plot 4 selected threshold values
ggplot(gcPlot, aes(x=Diversity, y=funcMaxed))+
  geom_point()+
  facet_wrap(~percent)+
  stat_smooth(method="glm", family = quasipoisson(link="identity"),colour="red", lwd=1.2 )+
  labs(y = expression("Number of Functions" >= Threshold), x = ("Species Richness"))+
  theme_bw(base_size=15)

# plot N Function over Threshold ~ Richness Slopes for all Thresholds
germanyThresh$percent <- 100*germanyThresh$thresholds 

ggplot(data=germanyThresh, aes(x=Diversity, y=funcMaxed, group=percent)) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("Species Richness") +
  stat_smooth(method="glm", family=quasipoisson(link="identity"), lwd=0.8,
              fill=NA, aes(color=percent)) +
  theme_bw(base_size=14) +
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red")

# plot the slopes of the relationship against Threshold values
germanyLinearSlopes<-getCoefTab(funcMaxed ~ Diversity, data=germanyThresh, coefVar="Diversity", family=quasipoisson(link="identity"))

ggplot(germanyLinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate- 1.96*germanyLinearSlopes[["Std. Error"]],
                                 ymax=Estimate+1.96*germanyLinearSlopes[["Std. Error"]])) + 
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Change in Number of Functions per Addition of 1 Species\n") + xlab("\nThreshold (%)") +
  stat_abline(intercept=0, slope=0, lwd=1, linetype=2) + 
  theme_bw(base_size=14)

# highlight critical slopes
germanIDX <- getIndices(germanyLinearSlopes, germanyThresh, funcMaxed ~ Diversity)
germanyThresh$IDX <- 0
germanyThresh$IDX [which(germanyThresh$thresholds %in% c(germanIDX$Tmin, germanIDX$Tmax,
                           germanIDX$Tmde))] <- 1

ggplot(germanyThresh, aes(x=Diversity, y=funcMaxed, group=percent)) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("Species Richness") +
  geom_smooth(method="glm", family=quasipoisson(link="identity"),
                fill=NA, aes(color=percent, lwd=IDX)) + theme_bw(base_size=14) +
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red") + 
  scale_size(range=c(0.3,5), guide="none") +
  annotate(geom="text", x=0, y=c(0.2,2,4.6), label=c("Tmax", "Tmde", "Tmin"))+
  annotate(geom="text", x=16.7, y=c(germanIDX$Mmin, germanIDX$Mmax,germanIDX$Mmde),
           label=c("Mmin", "Mmax", "Mmde"))

