###set your working dir (the badday2)

##### setting the right working dir
getwd()
bad2dir='/Users/luisacutillo/Sheffield/BAD DAYS/BAD_DAY2/'
setwd(bad2dir)
getwd()

#data=read.delim('./Data/growth.txt',header=TRUE)
data1=read.csv('./Data/growth1.csv')

dim(data1)

colnames(data1)

summary(data1)


attach(data1)


Treatment

data.aov = aov(Growth ~ Treatment)

#see the results: not very clear
data.aov

#We need to view the summary to have more information:

summary(data.aov)

TukeyHSD(data.aov)

data2=read.csv('./Data/growth2.csv')
names(data2)
detach(data1)




attach(data2)

Age

Treatment

length(Age)

tapply(Growth,Treatment,mean)
boxplot(Growth ~ Treatment)
tapply(Growth,Age,mean)
boxplot(Growth ~ Age)


int <- aov(Growth ~ Age*Treatment)
summary(int)

noint <- aov(Growth ~ Age+Treatment)
summary(noint)

library(MASS)

npk
replications(yield ~ N*P*K, data=npk)

#tapply
#?npk

with(npk, tapply(yield, list(N,P), mean));
with(npk, tapply(yield, list(N,K), mean));
with(npk, tapply(yield, list(P,K), mean));

npk.aov <- aov(yield ~ N*P*K, data=npk);
TukeyHSD(npk.aov, conf.level=.99);
plot(TukeyHSD(npk.aov, conf.level=.99))

summary(npk.aov)

plot(npk.aov);
#plot.design(yield~N*P*K, data=npk);
#qqnorm(npk$yield); qqline(npk$yield, col=4)



plot.design(yield~N*P*K, data=npk);
qqnorm(npk$yield); qqline(npk$yield, col=4)

by(npk$yield, npk$N, shapiro.test);
by(npk$yield, npk$P, shapiro.test);
by(npk$yield, npk$K, shapiro.test)

bartlett.test(npk$yield ~ npk$N)
bartlett.test(npk$yield ~ npk$P)
bartlett.test(npk$yield ~ npk$K)

fligner.test(npk$yield ~ npk$N);
fligner.test(npk$yield ~ npk$P);
fligner.test(npk$yield ~ npk$K);


