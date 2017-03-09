install.packages("xlsx")
install.packages("gdata")
install.packages("ape")

source("http://bioconductor.org/biocLite.R")
biocLite("multtest")

#################Exploratory plots

require(multtest)
data(golub)
?golub

#Gene expression dataset from Golub et al. (1999)
#Description

#Gene expression data (3051 genes and 38 tumor mRNA samples) from the leukemia 
#microarray study of Golub et al. (1999). Pre-processing was done as described 
#in Dudoit et al. (2002). 
dim(golub)
head(golub)
#add 3051 gene names
row.names(golub)=golub.gnames[,3]
#numeric vector indicating the tumor class, 27 acute lymphoblastic leukemia (ALL) cases (code 0) 
#and 11 acute myeloid leukemia (AML) cases (code 1).
colnames(golub)=golub.cl
head(golub)

###SCATTER PLOTS
#plot the value of gene CCND3 in all 38 mRNA samples (M92287_at)-- SCATTER PLOT
mygene <- golub[1042,]
plot(mygene)
#better plot custom
plot(mygene, pch=15, col="red");

#gene expression plot between patient 1 (ALL) and patient 38 (AML) and add the diagonal lines to the plot: abline(a= , b= )

plot(golub[,1], golub[,38])
abline(a=0, b=1, col="red", lwd=3)
###Scatterplots to detect indipendence
new.frame()
mysamplist <- golub[,c(1:5)]
colnames(mysamplist)=c(1:5)
plot(as.data.frame(mysamplist),pch='.')

#Bar plot of 4 cycline genes expression values in 3 ALL and AML patients respectively
mygenelist <- golub[c(85,1042,1212,2240),c(1:3, 36:38)]
#Make a bar plot with basic graph barplot function: barplot( )
barplot(mygenelist)

#bar plot arguments:
#beside=: to be portrayed as stacked bars or juxtaposed bars
#horiz=: bars to be drawn vertically or horizontally
#ylim=: limits for the y axis
#xlim=: limits for the x axis
#col=: color choices
barplot(mygenelist,horiz=TRUE,col=(c(1:4)),legend=TRUE)

###you can also use barplot to represent the mean and the standard error
ALLmean <- rowMeans(golub[c(85,1042,1212,2240),c(1:27)])
AMLmean <- rowMeans(golub[c(85,1042,1212,2240),c(28:38)])
dataheight <- cbind(ALLmean, AMLmean)
barx <- barplot(dataheight, beside=T, horiz=F, col=c(1:4), ylim=c(-2,2.5))
ALLsd <- apply(golub[c(85,1042,1212,2240),c(1:27)], 1, sd)
nALL=length(c(1:27))
AMLsd <- apply(golub[c(85,1042,1212,2240),c(28:38)], 1, sd)
nAML=length(c(28:38))

datasd <- cbind(ALLsd, AMLsd)
datase <- cbind(ALLsd/sqrt(nALL), AMLsd/sqrt(nAML))


par(mfrow=c(1,2))

###plot with the sd
datasdend<-abs(dataheight)+abs(datasd)
datasdend[c(3,4),]=-datasdend[c(3,4),]
barx <- barplot(dataheight, beside=T, horiz=F, col=c(1:4), ylim=c(-2,2.5))

arrows(barx, dataheight, barx, datasdend, angle=90, lwd=5, length=0.15, col=c(1:4))

###plot with the se THE ERROR ASSOCIATE TO THE MEAN!!!
datasdend<-abs(dataheight)+abs(datase)
datasdend[c(3,4),]=-datasdend[c(3,4),]
barx <- barplot(dataheight, beside=T, horiz=F, col=c(1:4), ylim=c(-2,2.5))

arrows(barx, dataheight, barx, datasdend, angle=90, lwd=2, length=0.2, col=1)





###Frequency table
#Discrete data occur when the values naturally fall into categories. A fre-
#  quency table simply gives the number of occurrences within a category.
#Example 1.
#A gene consists of a sequence of nucleotides (A;C;G;T)

#The number of each nucleotide can be displayed in a frequency table or a Pie plot
#The code below illustrates how to read the sequence "X94991.1" of the species homo
#sapiens from GenBank, , to construct a pie from a frequency table of the four
#nucleotides.

library('ape')

v=read.GenBank(c("X94991.1"),as.character = TRUE)
pie(table(v$X94991.1))

###stripcharts
#An elementary method to visualize data is by using a so-called stripchart,
#by which the values of the data are represented as e.g. small boxes
#it is useful in combination with a factor that distinguishes members from
#different experimental conditions or patients groups

#data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
stripchart(golub[1042,] ~ gol.fac, method="jitter")

#Histograms
#Another method to visualize data is by dividing the range of data values into
#a number of intervals and to plot the frequency per interval as a bar. Such
#a plot is called a histogram.

par(mfrow=c(2,2))
hist(golub[1042, gol.fac=="ALL"])
hist(golub,breaks=10)
hist(golub[, gol.fac=="AML"],breaks=10)
hist(golub[, gol.fac=="ALL"],breaks=10)

###BOXPLOTS
#A popular method to display data is by
#drawing a box around the Ørst and the third quartile (a bold line segment                                                     for the median), and the smaller line segments (whiskers) for the smallest and
#the largest data values. Such a data display is known as a box-and-whisker
#plot


# Sort the values of one gene
x <- sort(golub[1042, gol.fac=="ALL"], decreasing = FALSE)
x[1:5]

#A view on the distribution of the expression values of the
#ALL and the AML patients on gene CCND3 Cyclin D3 can be obtained by
#constructing two separate boxplots adjacent to one another. 
par(mfrow=c(1,1))
boxplot(golub[1042,] ~ gol.fac)
hist(golub[1042,], col="green", border=4, breaks=6, freq=F)
lines(density(golub[1042,]), col="red", lwd=5)

#Observe the distribution of all gene expression values in all 38 patients

boxplot(golub, col=2, lwd=2, border="green", pch=18)

#To compute exact values for the quartiles we need a sequence running
#from 0.00 to 1.00 with steps equal to 0.25. 
pvec <- seq(0,1,0.25)
quantile(golub[1042, gol.fac=="ALL"],pvec)

#Outliers are data values laying far apart from the pattern set by the
#majority of the data values. The implementation in R of the 
#boxplot draws such outlier points separately as small circles. A data point
#x is defined (graphically and not statistically) as an outlier point if x < x0.25 -1.5(X0.75-x0.25) [x>x0.25 >1.5(X0.75-x0.25)]

###Q-Q Plots
#A method to visualize the distribution of gene expression values is by the
#so-called quantile-quantile (Q-Q) plot. In such a plot the quantiles of the
#gene expression values are displayed against the corresponding quantiles of
#the normal (bell-shaped).

qqnorm(golub[1042, gol.fac=="ALL"])
qqline(golub[1042, gol.fac=="ALL"])

######################IMPORT--EXPORT data
#Load tab-delimited data - "NeuralStemCellData.tab"


mydata<-read.delim( "./Data/NeuralStemCellData.tab.txt", row.names=1, header=T)

####try do do some exploratory analisis on this data!

# GvHD flow cytometry data
# Only extract the CD3 positive cells
mydata<-as.data.frame(mydata)
cor(mydata[,1],mydata[,2])
plot(mydata[,1],mydata[,3])
