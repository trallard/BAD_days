
# <font color="#4682B4"> BAD Day 1 : T-test </font>



## 1. Loading the data set


```R
require(multtest)
data(golub)
```

    Loading required package: multtest
    Warning message in library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE, :
    “there is no package called ‘multtest’”Warning message in data(golub):
    “data set ‘golub’ not found”

We will be usig the Gene Expression dataset from Golub et al (1999).

Gene expression data (3051 genes and 38 tumor mRNA samples) from the leukemia microarray study of Golub et al. (1999). Pre-processing was done as described in Dudoit et al. (2002). The R code for pre-processing is available in the file ../doc/golub.R.
Source: Golub et al. (1999). Molecular classification of cancer: class discovery and class prediction by gene expression monitoring, Science, Vol. 286:531-537. http://www-genome.wi.mit.edu/MPR/ .


```R
golub.expr <- golub

# preliminar view of the data
head(golub.expr)
```


    Error in eval(expr, envir, enclos): object 'golub' not found
    Traceback:



 `golub.names` is a matrix containing the names of the 3051 
 genes contained in `golub`. The three columns correspond to:
`index, ID and name`


```R
row.names(golub.expr) = golub.gnames[,3]
```

`golub.cl` is a numeric vector indicating the tumor class, 27 acute lymphoblastic leukemia (ALL) cases (code 0) and 11 acute myeloid leukemia (AML) cases (code 1).


```R
colnames(golub.expr) = golub.cl
```

Now we need to set the sample sizes


```R
n.ALL <- 27
n.AML <- 11

cancer.type <- c(rep('ALL', n.ALL), rep('AML', n.AML))
```

Adding the cancer type to the column name, for the display


```R
colnames(golub.expr) <- cancer.type
```


```R
head(golub.expr)
```

## 2. t.test with a single gene (step-by-step)


```R
g <- 347
```

Alternatively, you can select a gene randomly


```R
# g <- sample(1:nrow(golub.expr),1)
g.profile <- as.vector(as.matrix(golub.expr[g,]))
```

Draw a barplot with color-coded cancer type


```R
plot.col <- c('ALL'='lightblue', 'AML'='mediumpurple')

barplot(g.profile, main=paste("Golub (1999), gene", g), 
        col=plot.col[cancer.type])
legend('topright', c("ALL","AML"),col=plot.col[c("ALL","AML")],
       pch=15, bty="o", bg='white')
box()
```

Separating the data in two vectors


```R
sample.ALL <- g.profile[cancer.type=="ALL"]
sample.AML <- g.profile[cancer.type=="AML"]
```

**Compute manually ** the t test parameters (not necessary, just to practice!)

Estimate the population means


```R
mean.est.ALL <- mean(sample.ALL)
mean.est.AML <- mean(sample.AML)
```

Compute the sample standard deviations 
The sd() function automatically computes the estimate corrected with sqrt(n-1)


```R
sample.sd.ALL <- sd(sample.ALL) * sqrt((n.ALL-1)/n.ALL)
sample.sd.AML <- sd(sample.AML) * sqrt((n.AML-1)/n.AML)
```

Estimate the population standard deviation



```R
sd.est.ALL <- sd(sample.ALL)
sd.est.AML <- sd(sample.AML)
```

Estimate the standard deviations on the means


```R
sd.err.est.ALL <- sd(sample.ALL) / sqrt(n.ALL)
sd.err.est.AML <- sd(sample.AML) / sqrt(n.AML)
```

Estimate the standard deviation of the difference between two means, according to Student's formula


```R
diff.sd.est <- sqrt((n.ALL*sample.sd.ALL^2 + n.AML*sample.sd.AML^2) * (1/n.ALL + 1/n.AML) /(n.ALL+n.AML-2))
```

Compute t.obs


```R
d <- abs(mean.est.ALL - mean.est.AML)
t.obs.Student <- d / diff.sd.est
```

Compute the P- vale.
Since we are performing the two-tail test, the single-tail probability has to be multiplied by 3 in order to obtai the alpha risk


```R
P.val.Student <- 2 * pt(q = t.obs.Student, df = n.ALL + n.AML-2, lower.tail = F)
```

## 3. T-test the fast way

### This is what you should be doing... 

### 3.1 Apply the Student-Fischer t-test (this assumes that the two populations have equal variance).


```R
t.student <- t.test(sample.ALL,sample.AML, var.equal=TRUE)
print(t.student)
```

### 3.2 Apply the Welch t-test (this does not assume that the two populations have equal variance)


```R
t.welch <- t.test(sample.ALL,sample.AML, var.equal=FALSE)
print(t.welch) 
```
