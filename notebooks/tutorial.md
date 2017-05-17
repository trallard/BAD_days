---
layout: default
title: "Tutorial"
tags:
    - notebook

permalink: "Tutorial.html"
---
# <font color="#4682B4"> BAD Day 1: Tutorial  </font>

# 0. Source/install the needed packages


<br>
<font color ='#00bcd4'> In [1]: </font>

{% highlight R %}
# In case you need to install the packages
install.packages("xlsx")
install.packages("gdata")
install.packages("ape")

{% endhighlight %}

    
    The downloaded binary packages are in
    	/var/folders/1q/xdx6qpy905dbx01t7cfv36280000gn/T//RtmpePD3cW/downloaded_packages
    
    The downloaded binary packages are in
    	/var/folders/1q/xdx6qpy905dbx01t7cfv36280000gn/T//RtmpePD3cW/downloaded_packages
    
    The downloaded binary packages are in
    	/var/folders/1q/xdx6qpy905dbx01t7cfv36280000gn/T//RtmpePD3cW/downloaded_packages


<br>
<font color ='#00bcd4'> In [2]: </font>

{% highlight R %}
source("http://bioconductor.org/biocLite.R");
biocLite("multtest");
{% endhighlight %}

    Bioconductor version 3.4 (BiocInstaller 1.24.0), ?biocLite for help
    A new version of Bioconductor is available after installing the most recent
      version of R; see http://bioconductor.org/install
    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.2 (2016-10-31).
    Installing package(s) ‘multtest’


    
    The downloaded binary packages are in
    	/var/folders/1q/xdx6qpy905dbx01t7cfv36280000gn/T//RtmpePD3cW/downloaded_packages


    Old packages: 'devtools', 'R6'



# 1. Exploratory data analysis

We will be usig the Gene Expression dataset from **Golub et al (1999)**. The
gene expression data collected by Golub et al. (1999) are among the most
classical in bioinformatics. A selection of the set is called `golub` which is
contained in the `multtest` package loaded before.


The data consist of gene expression values of 3051 genes (rows) from 38 leukemia
patients Pre-processing was done as described in Dudoit et al. (2002). The R
code for pre-processing is available in the file ../doc/golub.R.

**Source**:
Golub et al. (1999). Molecular classification of cancer: class discovery and
class prediction by gene expression monitoring, Science, Vol. 286:531-537.
http://www-genome.wi.mit.edu/MPR/ .

<br>
<font color ='#00bcd4'> In [4]: </font>

{% highlight R %}
require(multtest);

# Usage
data(golub);

# If you need more information on the data set just
# uncomment the line below
# ?golub
{% endhighlight %}

Data set values:
- `golub`: matrix of gene expression levels for the 38 tumor mRNA samples, rows
correspond to genes (3051 genes) and columns to mRNA samples.
- `golub.cl`: numeric vector indicating the tumor class, 27 acute lymphoblastic
leukemia (ALL) cases (code 0) and 11 acute myeloid leukemia (AML) cases (code
1).
- `golub.names`: a matrix containing the names of the 3051 genes for the
expression matrix golub. The three columns correspond to the gene index, ID, and
Name, respectively.

<br>
<font color ='#00bcd4'> In [5]: </font>

{% highlight R %}
# Checking the dimension of the data
dim(golub)
{% endhighlight %}





<br>
<font color ='#00bcd4'> In [6]: </font>

{% highlight R %}
# we will have a look at the first rows contained in the data set
head(golub)
{% endhighlight %}





The gene names are collected in the matrix `golub.gnames` of which the columns
correspond to the gene index, ID, and Name, respectively.

<br>
<font color ='#00bcd4'> In [7]: </font>

{% highlight R %}
# Adding 3051 gene names
row.names(golub) = golub.gnames[,3]

head(golub)
{% endhighlight %}





Twenty seven patients are diagnosed as acute lymphoblastic leukemia (ALL) and
eleven as acute myeloid leukemia (AML). The tumor class is given by the numeric
vector golub.cl, where ALL is indicated by 0 and AML by 1.

<br>
<font color ='#00bcd4'> In [8]: </font>

{% highlight R %}
colnames(golub) = golub.cl

head(golub)
{% endhighlight %}





## Creating the exploratory plots

### 1.1\. Plotting the value of gene (CCND3) in all nRNA samples (M92287_at)

We shall first have a look at the expression values of a gener with manufacurer
name `M92278_at`, which is known in biology as "CCND3 Cyclin D3".

The expression values of this gene are collected in row 1042 of golub. To load
the data and to obtain the relevant information from row 1042 of golub.gnames,
use the following:

<br>
<font color ='#00bcd4'> In [9]: </font>

{% highlight R %}
mygene <- golub[1042, ]
{% endhighlight %}

The data has now been stored in the `golub` matrix. We will now plot the
expression values od the gene CCND3 Cyclin D3.

<br>
<font color ='#00bcd4'> In [10]: </font>

{% highlight R %}
plot(mygene)
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_16_0.png)


In the previous plot we just used the default plotting preferences within R base
plotting.We can do some improvements so that the plot is easily understood.

<br>
<font color ='#00bcd4'> In [11]: </font>

{% highlight R %}
plot(mygene, pch = 15, col = 'slateblue', ylab = 'Expression value of gene: CCND3')
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_18_0.png)


In this plot the vertical axis corresponds to the size of the expression values
and the horizontal axis the index of the patients.

### 1.2\. Gene expression between patient 1 (ALL) and patient 38 (AML)

<br>
<font color ='#00bcd4'> In [12]: </font>

{% highlight R %}
plot(golub[,1], golub[,38])
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_21_0.png)


Adding diagonal lines to the plot and changing axes labels


<br>
<font color ='#00bcd4'> In [13]: </font>

{% highlight R %}
plot(golub[,1], golub[,38], xlab = 'Patient 1 (ALL)', ylab = 'Patient 38 (AML)') 
abline(a = 0, b = 1, col = 'mediumpurple', lwd =3)
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_23_0.png)


### 1.3\. Scatter plots to detect independence


<br>
<font color ='#00bcd4'> In [14]: </font>

{% highlight R %}
mysamplist <- golub[, c(1:15)]
colnames(mysamplist) = c(1:15)
{% endhighlight %}

<br>
<font color ='#00bcd4'> In [15]: </font>

{% highlight R %}
plot(as.data.frame(mysamplist), pch='.')
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_26_0.png)


### 1.4\. Bar plot of 4 cyclin genes expression values in 3 ALL and AML patients

We will analyse the expression values of the `D13639_at, M92287_at, U11791_at,
Z36714_AT` genes in three randomly chosen AML and ALL patients

<br>
<font color ='#00bcd4'> In [27]: </font>

{% highlight R %}
mygenelist <- golub[c(85, 1042, 1212, 2240), c(1:3, 36:38)]

# having a look at the data set chosen
mygenelist
{% endhighlight %}





<br>
<font color ='#00bcd4'> In [28]: </font>

{% highlight R %}
barplot(mygenelist, legend = T)
box()
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_29_0.png)


In this case the patients are indicated on the `X` axis (0 and 1 respectively)
while the gene expression level is indicate on the `Y` axis.

We can make some improvements to the plots.
Let's have a look at the `barplot` arguments:

<br>
<font color ='#00bcd4'> In [19]: </font>

{% highlight R %}
?barplot
{% endhighlight %}

We are going to focus on only a few of the histgram arguments:
- `beside`: `TRUE` for the bars to be displayed as justapoxed bars, `FALSE` for
stacked bars
- `horiz` : `FALSE` bars displayed vertically with the first bar to the left,
`TRUE` bars are displayed horizontally with the first at the bottom.
- `ylim`, `xlim` :  limits for the y and x axes
- `col`: colour choices

<br>
<font color ='#00bcd4'> In [30]: </font>

{% highlight R %}
# custom colours 
colours = c('lightblue',  'mediumpurple', 'lavender', 'midnightblue')

barplot(mygenelist, horiz = TRUE, col = colours, legend = TRUE,
       ylab = 'Patient', xlab = 'Gene expression level')
box()
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_33_0.png)


In the plot above we presented the barplots horizontally and added some colours,
which makes it easier to understand the data presented.
You can also use the barplots to represent the mean and standard error which we
will be doing in the following sections.

### 1.5\. Plotting the mean

In the following we will compute the mean for the expression values of both the
ALL and AML patients. We will be using the same 4 cycline genes used in the
example above.

<br>
<font color ='#00bcd4'> In [32]: </font>

{% highlight R %}
# Calculating the mean of the chosen genes from patient 1 to 27 and 28 to 38
ALLmean <- rowMeans(golub[c(85,1042,1212,2240),c(1:27)])
AMLmean <- rowMeans(golub[c(85,1042,1212,2240),c(28:38)])

# Combining the mean matrices previously calculated
dataheight <- cbind(ALLmean, AMLmean)

# Plotting 
barx <- barplot(dataheight, beside=T, horiz=F, col= colours, ylim=c(-2,2.5), 
                legend = TRUE, ylab = 'Gene expression level')
box()
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/tutorial_files/tutorial_36_0.png)


### 1.6\. Adding error bars to the previous plot

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
# Calculating the SD
ALLsd <- apply(golub[c(85,1042,1212,2240),c(1:27)], 1, sd)
nALL=length(c(1:27))
AMLsd <- apply(golub[c(85,1042,1212,2240),c(28:38)], 1, sd)
nAML=length(c(28:38))

datasd <- cbind(ALLsd, AMLsd)
datase <- cbind(ALLsd/sqrt(nALL), AMLsd/sqrt(nAML))

{% endhighlight %}

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
# creating a panel of 2 plots displayed in 1 row
par(mfrow = c(1,2))

# Plot with the SD
datasdend<-abs(dataheight) + abs(datasd)
datasdend[c(3,4),] = - datasdend[c(3,4),]
barx <- barplot(dataheight, beside=T, horiz=F, col= colours, ylim=c(-2,2.5),
               main = 'Data +  SD')

arrows(barx, dataheight, barx, datasdend, angle=90, lwd=5, length=0.15, 
       col = colours)
box()

# Plot with the se: error associated to the mean!
datasdend<-abs(dataheight) + abs(datase)
datasdend[c(3,4),] = -datasdend[c(3,4),]
barx <- barplot(dataheight, beside=T, horiz=F, col =colours, ylim=c(-2,2.5),
               main = 'Data + SE')

arrows(barx, dataheight, barx, datasdend, angle=90, lwd=2, length=0.2,
       col = 1)
box()
{% endhighlight %}

## 2. Exploring some types of plots


### 2.1 Frequency table
Discrete data occur when the values naturally fall into categories. A frequency
table simply gives the number of occurrences within a category.

A gene consists of a sequence of nucleotides (A; C; G; T)

The number of each nucleotide can be displayed in a frequency table.

This will be illustrated by the Zyxin gene which plays an important role in cell
adhesion The accession number (X94991.1) of one of its variants can be found in
a data base like NCBI (UniGene). The code below illustrates how to read the
sequence ”X94991.1” of the species homo sapiens from GenBank, to construct a
pie from a frequency table of the four nucleotides .

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
library('ape')
{% endhighlight %}

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
v = read.GenBank(c("X94991.1"),as.character = TRUE)

pie(table(v$X94991.1), col = colours)
{% endhighlight %}

### 2.2 Stripcharts

An elementary method to visualize data is by using a so-called stripchart,
by which the values of the data are represented as e.g. small boxes
it is useful in combination with a factor that distinguishes members from
different experimental conditions or patients groups.

Once again we use the CCND3 Cyclin D3 data to generate the plots.

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
# data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

stripchart(golub[1042,] ~ gol.fac, method = "jitter", 
           col = c('slateblue', 'darkgrey'), pch = 16)

{% endhighlight %}

To produce two adjacent stripcharts one for the ALL and one for the AML
patients, we use the factor called `gol.fac`. From the above figure, it can be
observed that the CCND3 Cyclin D3 expression values of the ALL patients tend to
have larger expression values than those of the AML patient.

### 2.3 Histograms

Another method to visualize data is by dividing the range of data values into
a number of intervals and to plot the frequency per interval as a bar. Such
a plot is called a histogram.

We will now generate a histogram of the expression values of gene CCND

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
par(mfrow=c(2,2))

hist(golub[1042, gol.fac == "ALL"], 
     col = 'lightblue', border = 'white',
    main = 'Golub[1042], ALL', xlab = 'ALL')
box()

hist(golub,breaks = 10, 
    col = 'slateblue', border = 'white',
    main =  'Golub')
box()

hist(golub[, gol.fac == "AML"],breaks = 10, 
     col = 'slateblue', border = 'white',
    main = 'Golub, AML', xlab = 'AML')
box()

hist(golub[, gol.fac == "ALL"],breaks = 10,
     col = 'slateblue', border = 'white',
    main = 'Golub, ALL', xlab = 'ALL')
box()
{% endhighlight %}

### 2.3 Boxplots

A popular method to display data is by
drawing a box around the Ørst and the third quartile (a bold line segment
for the median), and the smaller line segments (whiskers) for the smallest and
the largest data values. Such a data display is known as a box-and-whisker
plot

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
# Sort the values of one gene
x <- sort(golub[1042, gol.fac=="ALL"], decreasing = FALSE)
x[1:5]
{% endhighlight %}

A view on the distribution of the gene expression values of the `ALL` and `AML`
patients on gene CCND3 Cyclin D3 can be obtained by construction two separate
boxplots adjacent to each other:

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
# Even though we are creating two boxplots we only need one major graph
par(mfrow=c(1,1))
boxplot(golub[1042,] ~ gol.fac, col = c('lightblue', 'mediumpurple'))

{% endhighlight %}

It can be observed that the gene expression values for ALL are larger than those
for AML. Furthermore, since the two sub-boxes around the median are more or less
equally wide, the data are quite symmetrically distributed around the median.

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}

hist(golub[1042,], col= 'lightblue', border= 'black', breaks= 6, freq= F)
lines(density(golub[1042,]), col= 'slateblue', lwd = 3)
box()
{% endhighlight %}

Now we can observe the distribution of all gene expressions values in all 38
patients

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
boxplot(golub, col= 'lightblue', lwd=2, border="black", pch=18)
{% endhighlight %}

To compute exact values for the quartiles we need a sequence running from 0 to 1
with increments in steps of 0.25

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
pvec <- seq(0, 1, 0.25)
quantile(golub[1042, gol.fac=='ALL'], pvec)
{% endhighlight %}

Outliers are data points lying far apart from the pattern set by the majority of
the data values. The implementation in R of the boxplot draws such outliers as
smalle circles.

A data point `x` is defined (graphically, not statistically) as an outlier point
if $$x < 0.25 x -1.5\left(0.75 x -0.25 x\right) [x>0.25x >1.5(0.75x-0.25x)]$$


### 2.4 Q-Q plots

A method to visualize the distribution of gene expression values is y the so-
called quantile-quantile (Q-Q) plots. In such a plot the quantiles of the gene
expression values are displayed against the corresponding quantiles of the
normal distribution (bell-shaped).


To produce a Q-Q plot of the ALL gene expression values of CCND3 Cyclin D3 one
may use the following.

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
qqnorm(golub[1042, gol.fac == 'ALL'])
qqline(golub[1042, gol.fac == 'ALL'], col = 'mediumpurple', lwd = 3)
{% endhighlight %}

It can be seen that most of the data points are on or near the straight line,
while a few others are further away. The above example illustrates a case where
the degree of non-normality is moderate so that a clear conclusion cannot be
drawn.


## 3. Loading tab-delimited data

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
mydata<-read.delim("./NeuralStemCellData.tab.txt", row.names=1, header=T)
{% endhighlight %}

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
class(mydata)
{% endhighlight %}

### Now try and do some exploratory analysis of your own on this data!


GvHD flow cytometry data

Only exract the CD3 positive cells


<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}
cor(mydata[,1],mydata[,2])
plot(mydata[,1],mydata[,3])
{% endhighlight %}

<br>
<font color ='#00bcd4'> In [None]: </font>

{% highlight R %}

{% endhighlight %}
