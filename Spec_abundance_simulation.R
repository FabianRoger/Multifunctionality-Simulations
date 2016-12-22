
SPECNAM <- vector(length = 10000L)
SPECNAM <- vapply(SPECNAM, function(x = "a") x = paste(sample(LETTERS, 10, replace = TRUE), collapse = ""), c("a"))

SpecMat <- data.frame(SPECNAM = SPECNAM, V1 = round(rlnorm(10000, log(8), 2)))
SpecMat <- SpecMat[SpecMat$V1 != 0, ]


Species <- rep(SpecMat$SPECNAM, SpecMat$V1)
Species_sample <- data.frame(table(sample(Species, round(1 * length(Species)))))

diversity(Species_sample$Freq)

ggplot(SpecMat, aes(x = V1))+
   geom_histogram()+
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

ggplot(Species_sample, aes(x = Freq))+
  geom_histogram()+
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))



mean(log(SpecMat$V1))

exp(diversity(SpecMat$V1))

DF <- data.frame( proportion = 100*seq(0.01,1,0.01), 
                  effN = 1,
                  richness = 1,
                  Chao = 1)

for (i in seq(0.01,1,0.01)) {
  Species_sample <- data.frame(table(sample(Species, round(i * length(Species)), replace = T)))
  DF$effN[100*i] <- exp(diversity(Species_sample$Freq))
  DF$richness[100*i] <- specnumber(Species_sample$Freq)
  DF$Chao[100*i] <- ChaoSpecies(Species_sample$Freq)$Estimator
  
}

ggplot(DF, aes(x = proportion, y = Chao))+
  geom_point()+
  geom_point(data = DF, aes(x = proportion, y = richness), colour = "red")+
  


