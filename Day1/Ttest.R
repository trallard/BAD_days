####LOAD THE DATASET
require(multtest)
data(golub)


#Gene expression data (3051 genes and 38 tumor mRNA samples) from the leukemia 
#microarray study of Golub et al. (1999). Pre-processing was done as described 
#in Dudoit et al. (2002). 
golub.expr<-golub
row.names(golub.expr)=golub.gnames[,3]
#numeric vector indicating the tumor class, 27 acute lymphoblastic leukemia (ALL) cases (code 0) 
#and 11 acute myeloid leukemia (AML) cases (code 1).
colnames(golub.expr)=golub.cl

#####we need to set the nample sizes
n.ALL <- 27
n.AML <- 11
cancer.type <- c(rep("ALL", n.ALL), rep("AML", n.AML))

## Add the cancer type to the column name, for the display
colnames(golub.expr) <- cancer.type
# ##############################################################
## t.test with a single gene

g <- 347 

## Alternatively, you can select a gene randomly
## g <- sample(1:nrow(golub.expr),1)
g.profile <- as.vector(as.matrix(golub.expr[g,]))

## Draw a barpplot with color-coded cancer type
plot.col <- c('ALL'='#4444BB', 'AML'='#FFFF88')
x11(width=12,height=4)
barplot(g.profile,main=paste("Golub (1999), gene", g), col=plot.col[cancer.type])
legend('topright', c("ALL","AML"),col=plot.col[c("ALL","AML")],pch=15,bty="o",bg='white')
legend('topright', c("ALL","AML"),col='black',pch=22,bty="n")

## separate data in two vectors
sample.ALL <- g.profile[cancer.type=="ALL"]
sample.AML <- g.profile[cancer.type=="AML"]

#####COMPUTE MANUALLY the t test paramenters (not necessary, just to practice!)
## Estimate the population means
mean.est.ALL <- mean(sample.ALL)
mean.est.AML <- mean(sample.AML)

## Compute the sample sd
## Don't forget that the sd() function automatically computes the
## estimate corrected with sqrt(n-1)
sample.sd.ALL <- sd(sample.ALL) * sqrt((n.ALL-1)/n.ALL)
sample.sd.AML <- sd(sample.AML) * sqrt((n.AML-1)/n.AML)

## Estimate the population standard deviations 
## Don't forget that the sd() function automatically computes the
## estimate corrected with sqrt(n-1)
sd.est.ALL <- sd(sample.ALL)
sd.est.AML <- sd(sample.AML)

## Estimate the standard errors on the means
sd.err.est.ALL <- sd(sample.ALL) / sqrt(n.ALL)
sd.err.est.AML <- sd(sample.AML) / sqrt(n.AML)

## Estimate the standard deviation of the difference between two means, 
## according to Student's formula
diff.sd.est <- sqrt((n.ALL*sample.sd.ALL^2 + n.AML*sample.sd.AML^2) * (1/n.ALL + 1/n.AML) /(n.ALL+n.AML-2))

## Compute t.obs
d <- abs(mean.est.ALL - mean.est.AML)
t.obs.Student <- d/diff.sd.est

## Compute the P-value.
## Since we perform the two-tail test, the single-tail probability has
## to be multiplied by 2 in order to obtain the alpha risk.
P.val.Student <- 2*pt(q=t.obs.Student, df=n.ALL+n.AML-2,lower.tail=F)

#####THIS IS WHAT YOU SHOULD DO...FAST!
## Apply the Student-Fischer t-test 
## (this assumes that the two populations have an equal variance)
t.student <- t.test(sample.ALL,sample.AML, var.equal=TRUE)
print(t.student) 

## Apply the Welch t-test (
## (this does not assumes that the two populations 
## have an equal variance)
t.welch <- t.test(sample.ALL,sample.AML, var.equal=FALSE)
print(t.welch) 
