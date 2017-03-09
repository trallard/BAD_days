
#for more examples see: http://plantecology.syr.edu/fridley/bio793/glm.html
getwd()
###set your working dir (the badday2)

setwd('/Users/luisacutillo/Sheffield//BAD DAYS/BAD_DAY2/')
getwd()

dat = read.csv('./Data/treedata.csv') #choose the treedata.csv dataset
head(dat)
dim(dat)
dat2 = subset(dat,dat$species=="Tsuga canadensis")
mean(dat2$cover)



var(dat2$cover)


table(dat2$cover)


#If these counts were distributed exactly from a Poisson process, what would they look like, assuming the same mean (and variance)?

bar1 = barplot(as.vector(table(dat2$cover)),names.arg=seq(1:10))

points(bar1,dpois(seq(1,10),4.66)*sum(table(dat2$cover)),

     cex=2,type="b",col="sienna",lwd=2,pch=19)
#?dpois


#?glm

glm1 = glm(cover~1,data=dat2,family=poisson)

summary(glm1) 

with(dat2,plot(elev,cover,main="Hemlock cover vs. elevation",      cex=1.5,col="cyan",pch=19))
glm2 = glm(cover~elev,data=dat2,family=poisson)
summary(glm2)
anova(glm1,glm2,test="Chisq") 

#Note the chi-squared test is typically recommended for models with 'known deviance' (Poisson and binomial).  Here the model with elevation adds no explanatory power (fairly obvious from the graph), but we can still add the predicted trend line to our graph:

x = seq(0,1660)
#plot.new()
#lines(predict(glm2,list(elev=x),type="response"),lwd=2,col="orange")

#What is crucial here is the type argument to predict: "response" re-calculates the coefficients to be on the same scale as the original response variable, rather than the scale of the link function.  Let's now try an ANOVA with Poisson error, using disturbance as our predictor:

glm3 = glm(cover~disturb,data=dat2,family=poisson)

summary(glm3)
anova(glm1,glm3,test="Chisq") 


glm4 = glm(cover~disturb*elev,data=dat2,family=poisson)

summary(glm4) 

step(glm4) 

#?step



#st=step(glm4) 

#st$anova
