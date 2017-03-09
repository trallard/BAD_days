###show the application of the ANOVA (one way) on a built in dataset
help(PlantGrowth)
PlantGrowth
summary(PlantGrowth)
##compute sd of a given attribute
sd(PlantGrowth$weight) 
##compute sd of each group for the given attribute
tapply(PlantGrowth$weight, PlantGrowth$group, sd)
#boxplot of the data
boxplot(weight ~ group, data=PlantGrowth, col="blue")
#draw a conditional histogram
par(mfrow=c(1,3))
library(lattice) 
histogram(~ weight | group, data=PlantGrowth, col="blue", type="count") 

#Our null hypothesis is that the three groups have the same growth mean
analysis <- aov(weight ~ group, data=PlantGrowth)
summary(analysis) 
#build a linear model (using the function lm() command) and then use the anova() command to analyse that
fit <- lm(weight ~ group, data=PlantGrowth)
anova(fit) 
###same result!

