###set your working dir (the badday2)
bad2dir='/Users/luisacutillo/Sheffield/BAD DAYS/BAD_DAY2/'
setwd(bad2dir)
getwd()


install.packages(rafalib)
install.packages(RColorBrewer)
install.packages(Rtsne)
install.packages(tsne)
install.packages(pracma)
install.packages(Cairo)



data(iris)
str(iris); 
summary(iris[1:4])
#dim(iris)
####scatter plot of dimensions against each other
pairs(iris[1:4],main="Iris Data", pch=19, col=as.numeric(iris$Species)+1)

mtext("Color map for iris species: red-> setosa; green-> versicolor; blue-> virginica", 1, line=3.7,cex=.8)

#To examine variability of all numeric variables
sapply(iris[1:4],var)
range(sapply(iris[1:4],var))
# maybe this range of variability is big in this context.
#Thus, we will use the correlation matrix
#For this, we must standardize our variables with scale() function:
iris.stand <- as.data.frame(scale(iris[,1:4]))
sapply(iris.stand,sd) #now, standard deviations are 1
#sapply(iris.stand,mean) #now, standard deviations are 1


#If we use prcomp() function, we indicate 'scale=TRUE' to use correlation matrix
pca <- prcomp(iris.stand,scale=T)

#similar with princomp(): princomp(iris.stand, cor=T)
pca
summary(pca)
#This gives us the standard deviation of each component, and the proportion of variance explained by each component.
#The standard deviation is stored in (see 'str(pca)'):
pca$sdev

#plot of variance of each PCA.
#It will be useful to decide how many principal components should be retained.
screeplot(pca, type="lines",col=3)

#retreive loadings for the principal components:
pca$rotation # when using princomp(): pca$loadings

unique(iris$Species)
as.numeric(unique(iris$Species))+1

#biplot of first two principal components
biplot(pca,cex=0.8)
abline(h = 0, v = 0, lty = 2, col = 8)
plot(pca$x,col=as.numeric(iris$Species)+1); 
legend("topright",legend=unique(iris$Species),col=unique(as.numeric(iris$Species)+1),
         ,pch=1)

####fist of all, just as an example, apply tsne on the IRIS dataset
tsne_iris = tsne(iris[,1:4])

plot(tsne_iris,col=as.numeric(iris$Species)+1); 
legend("topright",legend=unique(iris$Species),col=unique(as.numeric(iris$Species)+1),
         ,pch=1)


#install flowCore package from Bioconductor (to read FCS files)
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")

# install Rtsne package from CRAN (R implementation of Barnes-Hut-SNE algorithm)

install.packages("Rtsne")

# load packages

library(flowCore)
library(Rtsne)


getwd()

data <- exprs(read.FCS("./Rtsne-example-master/data/visne_marrow1.fcs", transformation = FALSE))
head(data)

unname(colnames(data))  # isotope and marker (protein) names




# select markers to use in calculation of t-SNE projection
# CD11b, CD123, CD19, CD20, CD3, CD33, CD34, CD38, CD4, CD45, CD45RA, CD8, CD90
# (see Amir et al. 2013, Supplementary Tables 1 and 2)

colnames_proj <- unname(colnames(data))[c(11, 23, 10, 16, 7, 22, 14, 28, 12, 6, 8, 13, 30)]
colnames_proj  # check carefully!


# arcsinh transformation
# (see Amir et al. 2013, Online Methods, "Processing of mass cytometry data")

asinh_scale <- 5
data <- asinh(data / asinh_scale)  # transforms all columns! including event number etc


# subsampling

nsub <- 10000
set.seed(123)  # set random seed
data <- data[sample(1:nrow(data), nsub), ]

dim(data)


# prepare data for Rtsne

data <- data[, colnames_proj]      # select columns to use
data <- data[!duplicated(data), ]  # remove rows containing duplicate values within rounding

dim(data)



file <- paste("./Rtsne-example-master/data/viSNE_Marrow1_nsub", nsub, ".txt", sep = "")
write.table(data, file = file, row.names = FALSE, quote = FALSE, sep = "\t")


# run Rtsne (Barnes-Hut-SNE algorithm)
# without PCA step (see Amir et al. 2013, Online Methods, "viSNE analysis")

set.seed(123)  # set random seed
rtsne_out <- Rtsne(as.matrix(data), pca = FALSE, verbose = TRUE)


# plot 2D t-SNE projection

file_plot <- paste("./Rtsne-example-master/plots/Rtsne_viSNE_Marrow1_nsub", nsub, ".png", sep = "")
png(file_plot, width = 900, height = 900)
plot(rtsne_out$Y, asp = 1, pch = 20, col = "blue", 
     cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, 
     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", 
     main = "2D t-SNE projection")
dev.off()



plot(rtsne_out$Y, asp = 1, pch = 20, col = "blue", 
     cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, 
     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", 
     main = "2D t-SNE projection")








# compare to PCA
pca_data = princomp(as.matrix(data))$scores[,1:2]
plot(pca_data,col="blue")
#text(pca_iris, labels=iris$Species,col=colors[iris$Species])



library(rafalib)
library(RColorBrewer)
library(Rtsne)
library(tsne)
library(pracma)
library(Cairo)

par(mfrow=c(1, 2))
library(Rtsne)
library(tsne)
for (i in 2:5) {
  set.seed(i)
  x <- runif(n, -1, 1)
  cols <- brewer.pal(11, "Spectral")[as.integer(cut(x, 11))]
  ortho <- rortho(m)
  X <- cbind(x, matrix(0,ncol=m-1,nrow=n)) %*% ortho
  res <- tsne(X)
  plot(res, col=cols, pch=20, xlab="", ylab="", main="tsne")
  res <- Rtsne(X)
  plot(res$Y, col=cols, pch=20, xlab="", ylab="", main="Rtsne")
}

n <- 200
m <- 40
set.seed(1)
x <- runif(n, -1, 1)
library(rafalib)
bigpar(2,2,mar=c(3,3,3,1))
library(RColorBrewer)
cols <- brewer.pal(11, "Spectral")[as.integer(cut(x, 11))]
plot(x, rep(0,n), ylim=c(-1,1), yaxt="n", xlab="", ylab="",
col=cols, pch=20, main="underlying data")
library(pracma)
ortho <- rortho(m)
X <- cbind(x, matrix(0,ncol=m-1,nrow=n)) %*% ortho
plot(X[,1:2], asp=1, col=cols, pch=20, xlab="", ylab="", main="embed in higher dim")
pc <- prcomp(X)
plot(pc$x[,1:2], asp=1, col=cols, pch=20, xlab="", ylab="", main="PC1 & PC2")
library(tsne)
res <- tsne(X)
plot(res, col=cols, pch=20, xlab="", ylab="", main="t-SNE")
bigpar(2,2,mar=c(3,3,1,1))
for (i in 2:5) {
set.seed(i)
x <- runif(n, -1, 1)
cols <- brewer.pal(11, "Spectral")[as.integer(cut(x, 11))]
ortho <- rortho(m)
X <- cbind(x, matrix(0,ncol=m-1,nrow=n)) %*% ortho
res <- tsne(X)
plot(res, col=cols, pch=20, xlab="", ylab="")
}


