---
layout: default
title: "RNA-sequencing"

tags:
    - Day3
nb: "RNA-sequencing.ipynb"
module: '/3-outline/'

permalink: "RNA-sequencing.html"
---
# BAD Day 3: RNA sequencing pipeline

# BAD Day 3: RNA sequencing pipeline

This tutorial follows the Work flow described in
<http://www.bioconductor.org/help/workflows/RNAseq123/>.
Note that to install this workflow under Bioconductor 3.5 you need to run the
following commands:

```R
source("http://bioconductor.org/workflows.R")
workflowInstall("RNAseq123")
```


## Introduction

RNA-sequencing (RNA-seq) has become the primary technology used for gene
expression profiling, with the genome-wide detection of differentially expressed
genes between two or more conditions of interest one of the most commonly asked
questions by researchers. The edgeR (Robinson, McCarthy, and Smyth 2010) and
limma packages (Ritchie et al. 2015) available from the Bioconductor project
(Huber et al. 2015) offer a well-developed suite of statistical methods for
dealing with this question for RNA-seq data.

In this article, we describe an edgeR - limma workflow for analysing RNA-seq
data that takes gene-level counts as its input, and moves through pre-processing
and exploratory data analysis before obtaining lists of differentially expressed
(DE) genes and gene signatures. This analysis is enhanced through the use of
interactive graphics from the Glimma package (Su and Ritchie 2016), that allows
for a more detailed exploration of the data at both the sample and gene-level
than is possible using static R plots.

The experiment analysed in this workflow is from Sheridan et al. (2015)
(Sheridan et al. 2015) and consists of three cell populations (basal, luminal
progenitor (LP) and mature luminal (ML)) sorted from the mammary glands of
female virgin mice, each profiled in triplicate. RNA samples were sequenced
across three batches on an Illumina HiSeq 2000 to obtain 100 base-pair single-
end reads. The analysis outlined in this article assumes that reads obtained
from an RNA-seq experiment have been aligned to an appropriate reference genome
and summarised into counts associated with gene-specific regions. In this
instance, reads were aligned to the mouse reference genome (mm10) using the R
based pipeline available in the Rsubread package (specifically the align
function (Liao, Smyth, and Shi 2013) followed by featureCounts (Liao, Smyth, and
Shi 2014) for gene-level summarisation based on the in-built mm10 RefSeq-based
annotation).

Count data for these samples can be downloaded from the Gene Expression Omnibus
(GEO) http://www.ncbi.nlm.nih.gov/geo/ using GEO Series accession number
GSE63310. Further information on experimental design and sample preparation is
also available from GEO under this accession number.

<br>
<font color ='#00bcd4'> In [1]: </font>

{% highlight R %}
source("http://bioconductor.org/workflows.R")
workflowInstall("RNAseq123")

{% endhighlight %}

<br>
<font color ='#00bcd4'> In [2]: </font>

{% highlight R %}
install.packages("R.utils")

{% endhighlight %}

## Loading the required packages

<br>
<font color ='#00bcd4'> In [3]: </font>

{% highlight R %}
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus);
{% endhighlight %}

## Getting the data

The dataset ` GSE63310_RAW.tar ` will be directly donwloaded from
<https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file>. We will
then extract the relevant files.

<br>
<font color ='#00bcd4'> In [4]: </font>

{% highlight R %}
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"

utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb")
utils::untar("GSE63310_RAW.tar", exdir = ".")

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")

for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)
{% endhighlight %}

Each of these text files contains the raw gene-level counts for a given sample.

This analysis only includes the basal, LP and ML samples from this experiment
(see associated file names below).


<br>
<font color ='#00bcd4'> In [5]: </font>

{% highlight R %}
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",
   "GSM1545545_JMS9-P8c.txt")

read.delim(files[1], nrow=5)
files
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th scope="col">EntrezID</th><th scope="col">GeneLength</th><th scope="col">Count</th></tr></thead>
<tbody>
	<tr><td>   497097</td><td>3634     </td><td>1        </td></tr>
	<tr><td>100503874</td><td>3259     </td><td>0        </td></tr>
	<tr><td>100038431</td><td>1634     </td><td>0        </td></tr>
	<tr><td>    19888</td><td>9747     </td><td>0        </td></tr>
	<tr><td>    20671</td><td>3130     </td><td>1        </td></tr>
</tbody>
</table>




<ol class="list-inline">
	<li>'GSM1545535_10_6_5_11.txt'</li>
	<li>'GSM1545536_9_6_5_11.txt'</li>
	<li>'GSM1545538_purep53.txt'</li>
	<li>'GSM1545539_JMS8-2.txt'</li>
	<li>'GSM1545540_JMS8-3.txt'</li>
	<li>'GSM1545541_JMS8-4.txt'</li>
	<li>'GSM1545542_JMS8-5.txt'</li>
	<li>'GSM1545544_JMS9-P7c.txt'</li>
	<li>'GSM1545545_JMS9-P8c.txt'</li>
</ol>



Instead of readig each file separately and combine posteriorly the package
`edgeR` allows us to do this in one single step using the function `readDGE`.

The resulting object will be a matrix of counts with 27,179 rows associated with
unique Entrez gene identifiers (IDs) and nine columns associated with the
individual samples in the experiment.

<br>
<font color ='#00bcd4'> In [6]: </font>

{% highlight R %}
x <- readDGE(files, columns=c(1,3))
class(x)
{% endhighlight %}


'DGEList'


<br>
<font color ='#00bcd4'> In [7]: </font>

{% highlight R %}
# checking the dimension of x
dim(x)
x$samples
x$counts[1,]
{% endhighlight %}


<ol class="list-inline">
	<li>27179</li>
	<li>9</li>
</ol>




<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">files</th><th scope="col">group</th><th scope="col">lib.size</th><th scope="col">norm.factors</th></tr></thead>
<tbody>
	<tr><th scope="row">GSM1545535_10_6_5_11</th><td>GSM1545535_10_6_5_11.txt</td><td>1                       </td><td>32863052                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545536_9_6_5_11</th><td>GSM1545536_9_6_5_11.txt </td><td>1                       </td><td>35335491                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545538_purep53</th><td>GSM1545538_purep53.txt  </td><td>1                       </td><td>57160817                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545539_JMS8-2</th><td>GSM1545539_JMS8-2.txt   </td><td>1                       </td><td>51368625                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545540_JMS8-3</th><td>GSM1545540_JMS8-3.txt   </td><td>1                       </td><td>75795034                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545541_JMS8-4</th><td>GSM1545541_JMS8-4.txt   </td><td>1                       </td><td>60517657                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545542_JMS8-5</th><td>GSM1545542_JMS8-5.txt   </td><td>1                       </td><td>55086324                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545544_JMS9-P7c</th><td>GSM1545544_JMS9-P7c.txt </td><td>1                       </td><td>21311068                </td><td>1                       </td></tr>
	<tr><th scope="row">GSM1545545_JMS9-P8c</th><td>GSM1545545_JMS9-P8c.txt </td><td>1                       </td><td>19958838                </td><td>1                       </td></tr>
</tbody>
</table>




<dl class="dl-horizontal">
	<dt>GSM1545535_10_6_5_11</dt>
		<dd>1</dd>
	<dt>GSM1545536_9_6_5_11</dt>
		<dd>2</dd>
	<dt>GSM1545538_purep53</dt>
		<dd>342</dd>
	<dt>GSM1545539_JMS8-2</dt>
		<dd>526</dd>
	<dt>GSM1545540_JMS8-3</dt>
		<dd>3</dd>
	<dt>GSM1545541_JMS8-4</dt>
		<dd>3</dd>
	<dt>GSM1545542_JMS8-5</dt>
		<dd>535</dd>
	<dt>GSM1545544_JMS9-P7c</dt>
		<dd>2</dd>
	<dt>GSM1545545_JMS9-P8c</dt>
		<dd>0</dd>
</dl>



## Organising sample information


For downstream analysis, sample-level information related to the experimental
design needs to be associated with the columns of the counts matrix. This should
include experimental variables, both biological and technical, that could have
an effect on expression levels. Examples include cell type (basal, LP and ML in
this experiment), genotype (wild-type, knock-out), phenotype (disease status,
sex, age), sample treatment (drug, control) and batch information (date
experiment was performed if samples were collected and analysed at distinct time
points) to name just a few.

Our DGEList-object contains a `samples` data frame that stores both cell type
(or `group`) and batch (sequencing `lane`) information, each of which consists
of three distinct levels. Note that within `x$samples`, library sizes are
automatically calculated for each sample and normalisation factors are set to 1.
For simplicity, we remove the GEO sample IDs (GSM*) from the column names of our
DGEList-object `x`.

<br>
<font color ='#00bcd4'> In [8]: </font>

{% highlight R %}
colnames(x) #these are the file names
{% endhighlight %}


<ol class="list-inline">
	<li>'GSM1545535_10_6_5_11'</li>
	<li>'GSM1545536_9_6_5_11'</li>
	<li>'GSM1545538_purep53'</li>
	<li>'GSM1545539_JMS8-2'</li>
	<li>'GSM1545540_JMS8-3'</li>
	<li>'GSM1545541_JMS8-4'</li>
	<li>'GSM1545542_JMS8-5'</li>
	<li>'GSM1545544_JMS9-P7c'</li>
	<li>'GSM1545545_JMS9-P8c'</li>
</ol>



<br>
<font color ='#00bcd4'> In [9]: </font>

{% highlight R %}
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

{% endhighlight %}


<ol class="list-inline">
	<li>'10_6_5_11'</li>
	<li>'9_6_5_11'</li>
	<li>'purep53'</li>
	<li>'JMS8-2'</li>
	<li>'JMS8-3'</li>
	<li>'JMS8-4'</li>
	<li>'JMS8-5'</li>
	<li>'JMS9-P7c'</li>
	<li>'JMS9-P8c'</li>
</ol>



<br>
<font color ='#00bcd4'> In [10]: </font>

{% highlight R %}
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP",
                     "Basal", "ML", "LP"))

x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))

x$samples$lane <- lane
x$samples
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">files</th><th scope="col">group</th><th scope="col">lib.size</th><th scope="col">norm.factors</th><th scope="col">lane</th></tr></thead>
<tbody>
	<tr><th scope="row">10_6_5_11</th><td>GSM1545535_10_6_5_11.txt</td><td>LP                      </td><td>32863052                </td><td>1                       </td><td>L004                    </td></tr>
	<tr><th scope="row">9_6_5_11</th><td>GSM1545536_9_6_5_11.txt </td><td>ML                      </td><td>35335491                </td><td>1                       </td><td>L004                    </td></tr>
	<tr><th scope="row">purep53</th><td>GSM1545538_purep53.txt  </td><td>Basal                   </td><td>57160817                </td><td>1                       </td><td>L004                    </td></tr>
	<tr><th scope="row">JMS8-2</th><td>GSM1545539_JMS8-2.txt   </td><td>Basal                   </td><td>51368625                </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-3</th><td>GSM1545540_JMS8-3.txt   </td><td>ML                      </td><td>75795034                </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-4</th><td>GSM1545541_JMS8-4.txt   </td><td>LP                      </td><td>60517657                </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-5</th><td>GSM1545542_JMS8-5.txt   </td><td>Basal                   </td><td>55086324                </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS9-P7c</th><td>GSM1545544_JMS9-P7c.txt </td><td>ML                      </td><td>21311068                </td><td>1                       </td><td>L008                    </td></tr>
	<tr><th scope="row">JMS9-P8c</th><td>GSM1545545_JMS9-P8c.txt </td><td>LP                      </td><td>19958838                </td><td>1                       </td><td>L008                    </td></tr>
</tbody>
</table>



## Organising gene annotations


A second data frame named `genes` in the DGEList-object is used to store gene-
level information associated with rows of the counts matrix. This information
can be retrieved using organism specific packages such as `Mus.musculus`
(Bioconductor Core Team 2016b) for mouse (or `Homo.sapiens` (Bioconductor Core
Team 2016a) for human) or the `biomaRt` package (Durinck et al. 2005; Durinck et
al. 2009) which interfaces the Ensembl genome databases in order to perform gene
annotation.


The type of information that can be retrieved includes gene symbols, gene names,
chromosome names and locations, Entrez gene IDs, Refseq gene IDs and Ensembl
gene IDs to name just a few. `biomaRt` primarily works off Ensembl gene IDs,
whereas `Mus.musculus` packages information from various sources and allows
users to choose between many different gene IDs as the key.

The Entrez gene IDs available in our dataset were annotated using the
`Mus.musculus` package to retrieve associated gene symbols and chromosome
information.

<br>
<font color ='#00bcd4'> In [11]: </font>

{% highlight R %}
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"),
                keytype="ENTREZID")

head(genes)
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th scope="col">ENTREZID</th><th scope="col">SYMBOL</th><th scope="col">TXCHROM</th></tr></thead>
<tbody>
	<tr><td>497097   </td><td>Xkr4     </td><td>chr1     </td></tr>
	<tr><td>100503874</td><td>Gm19938  </td><td>NA       </td></tr>
	<tr><td>100038431</td><td>Gm10568  </td><td>NA       </td></tr>
	<tr><td>19888    </td><td>Rp1      </td><td>chr1     </td></tr>
	<tr><td>20671    </td><td>Sox17    </td><td>chr1     </td></tr>
	<tr><td>27395    </td><td>Mrpl15   </td><td>chr1     </td></tr>
</tbody>
</table>



<br>
<font color ='#00bcd4'> In [12]: </font>

{% highlight R %}
length(geneid)
length(genes$ENTREZID)
{% endhighlight %}


27179



27220


As with any gene ID, Entrez gene IDs may not map one-to-one to the gene
information of interest. It is important to check for duplicated gene IDs and to
understand the source of duplication before resolving them. Our gene annotation
contains 28 genes that map to multiple chromosomes (e.g. gene Gm1987 is
associated with chr4 and chr4_JH584294_random and microRNA Mir5098 is associated
with chr2, chr5, chr8, chr11 and chr17).

To resolve duplicate gene IDs one could combine all chromosome information from
the multi-mapped genes, such that gene Gm1987 would be is assigned to chr4 and
chr4_JH584294_random, or select one of the chromosomes to represent the gene
with duplicate annotation. For simplicity we do the latter, keeping only the
first occurrence of each gene ID.

<br>
<font color ='#00bcd4'> In [13]: </font>

{% highlight R %}
genes$ENTREZID[which(duplicated(genes$ENTREZID))]
{% endhighlight %}


<ol class="list-inline">
	<li>'100316809'</li>
	<li>'12228'</li>
	<li>'433182'</li>
	<li>'100217457'</li>
	<li>'735274'</li>
	<li>'735274'</li>
	<li>'735274'</li>
	<li>'735274'</li>
	<li>'735274'</li>
	<li>'735274'</li>
	<li>'100504346'</li>
	<li>'621580'</li>
	<li>'545611'</li>
	<li>'100040048'</li>
	<li>'100040048'</li>
	<li>'16158'</li>
	<li>'24047'</li>
	<li>'100316707'</li>
	<li>'100039939'</li>
	<li>'108816'</li>
	<li>'108816'</li>
	<li>'100504362'</li>
	<li>'100628577'</li>
	<li>'100628577'</li>
	<li>'100628577'</li>
	<li>'100628577'</li>
	<li>'381724'</li>
	<li>'331195'</li>
	<li>'331195'</li>
	<li>'100041102'</li>
	<li>'100041102'</li>
	<li>'622894'</li>
	<li>'100041253'</li>
	<li>'100041253'</li>
	<li>'100041354'</li>
	<li>'665943'</li>
	<li>'81016'</li>
	<li>'81017'</li>
	<li>'100039499'</li>
	<li>'30058'</li>
	<li>'100499528'</li>
</ol>



<br>
<font color ='#00bcd4'> In [14]: </font>

{% highlight R %}
genes$SYMBOL[which(duplicated(genes$ENTREZID))]
{% endhighlight %}


<ol class="list-inline">
	<li>'Mir1906-1'</li>
	<li>'Btg3'</li>
	<li>'Eno1b'</li>
	<li>'Snord58b'</li>
	<li>'Mir684-1'</li>
	<li>'Mir684-1'</li>
	<li>'Mir684-1'</li>
	<li>'Mir684-1'</li>
	<li>'Mir684-1'</li>
	<li>'Mir684-1'</li>
	<li>'Gm13304'</li>
	<li>'Gm13308'</li>
	<li>'Fam205a2'</li>
	<li>'Ccl27b'</li>
	<li>'Ccl27b'</li>
	<li>'Il11ra2'</li>
	<li>'Ccl19'</li>
	<li>'Mir1957a'</li>
	<li>'Gm2506'</li>
	<li>'4933409K07Rik'</li>
	<li>'4933409K07Rik'</li>
	<li>'Gm1987'</li>
	<li>'Mir5098'</li>
	<li>'Mir5098'</li>
	<li>'Mir5098'</li>
	<li>'Mir5098'</li>
	<li>'BC061212'</li>
	<li>'A430089I19Rik'</li>
	<li>'A430089I19Rik'</li>
	<li>'Gm3139'</li>
	<li>'Gm3139'</li>
	<li>'Gm6367'</li>
	<li>'Gm16513'</li>
	<li>'Gm16513'</li>
	<li>'Gm3286'</li>
	<li>'E330014E10Rik'</li>
	<li>'Vmn1r62'</li>
	<li>'Vmn1r63'</li>
	<li>'Vmn1r187'</li>
	<li>'Timm8a1'</li>
	<li>'Mir3473a'</li>
</ol>



<br>
<font color ='#00bcd4'> In [15]: </font>

{% highlight R %}
kk=c(1,1,1,2,2,3,4,4,5,5,5)
duplicated(kk)

{% endhighlight %}


<ol class="list-inline">
	<li>FALSE</li>
	<li>TRUE</li>
	<li>TRUE</li>
	<li>FALSE</li>
	<li>TRUE</li>
	<li>FALSE</li>
	<li>FALSE</li>
	<li>TRUE</li>
	<li>FALSE</li>
	<li>TRUE</li>
	<li>TRUE</li>
</ol>



<br>
<font color ='#00bcd4'> In [16]: </font>

{% highlight R %}
genes <- genes[!duplicated(genes$ENTREZID),]
{% endhighlight %}

In this example, the gene order is the same in both the annotation and the data
object. If this is not the case due to missing and/or rearranged gene IDs, the
`match` function can be used to order genes correctly. The data frame of gene
annotations is then added to the data object and neatly packaged in a DGEList-
object containing raw count data with associated sample information and gene
annotations.

<br>
<font color ='#00bcd4'> In [17]: </font>

{% highlight R %}
x$genes <- genes

{% endhighlight %}

## Data pre-processing

### Transformations from the raw-scale

For differential expression and related analyses, gene expression is rarely
considered at the level of raw counts since libraries sequenced at a greater
depth will result in higher counts. Rather, it is common practice to transform
raw counts onto a scale that accounts for such library size differences. Popular
transformations include counts per million (CPM), log2-counts per million (log-
CPM), reads per kilobase of transcript per million (RPKM), and fragments per
kilobase of transcript per million (FPKM).

In our analyses, CPM and log-CPM transformations are used regularly although
they do not account for feature length differences which RPKM and FPKM values
do. Whilst RPKM and FPKM values can just as well be used, CPM and log-CPM values
can be calculated using a counts matrix alone and will suffice for the type of
comparisons we are interested in. Assuming that there are no differences in
isoform usage between conditions, differential expression analyses look at gene
expression changes between conditions rather than comparing expression across
multiple genes or drawing conclusions on absolute levels of expression. In other
words, gene lengths remain constant for comparisons of interest and any observed
differences are a result of changes in condition rather than changes in gene
length.

Here raw counts are converted to CPM and log-CPM values using the `cpm` function
in `edgeR`, where log-transformations use a prior count of 0.25 to avoid taking
the log of zero. RPKM values are just as easily calculated as CPM values using
the `rpkm` function in `edgeR` if gene lengths are available.

<br>
<font color ='#00bcd4'> In [18]: </font>

{% highlight R %}
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

{% endhighlight %}

### Removing genes with low expression level

All datasets will include a mix of genes that are expressed and those that are
not expressed. Whilst it is of interest to examine genes that are expressed in
one condition but not in another, some genes are unexpressed throughout all
samples. In fact, 19% of genes in this dataset have zero counts across all nine
samples.

<br>
<font color ='#00bcd4'> In [19]: </font>

{% highlight R %}
table(rowSums(x$counts==0)==9)
keep.exprs <- rowSums(cpm>1)>=3

{% endhighlight %}



    FALSE  TRUE
    22026  5153


<br>
<font color ='#00bcd4'> In [20]: </font>

{% highlight R %}
x$samples$lib.size[1]
{% endhighlight %}


32863052


<br>
<font color ='#00bcd4'> In [21]: </font>

{% highlight R %}
x[c(1,2,3),,keep.lib.sizes=FALSE]
{% endhighlight %}


<dl>
	<dt>$samples</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">files</th><th scope="col">group</th><th scope="col">lib.size</th><th scope="col">norm.factors</th><th scope="col">lane</th></tr></thead>
<tbody>
	<tr><th scope="row">10_6_5_11</th><td>GSM1545535_10_6_5_11.txt</td><td>LP                      </td><td>  1                     </td><td>1                       </td><td>L004                    </td></tr>
	<tr><th scope="row">9_6_5_11</th><td>GSM1545536_9_6_5_11.txt </td><td>ML                      </td><td>  2                     </td><td>1                       </td><td>L004                    </td></tr>
	<tr><th scope="row">purep53</th><td>GSM1545538_purep53.txt  </td><td>Basal                   </td><td>347                     </td><td>1                       </td><td>L004                    </td></tr>
	<tr><th scope="row">JMS8-2</th><td>GSM1545539_JMS8-2.txt   </td><td>Basal                   </td><td>532                     </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-3</th><td>GSM1545540_JMS8-3.txt   </td><td>ML                      </td><td>  3                     </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-4</th><td>GSM1545541_JMS8-4.txt   </td><td>LP                      </td><td>  3                     </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-5</th><td>GSM1545542_JMS8-5.txt   </td><td>Basal                   </td><td>541                     </td><td>1                       </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS9-P7c</th><td>GSM1545544_JMS9-P7c.txt </td><td>ML                      </td><td>  2                     </td><td>1                       </td><td>L008                    </td></tr>
	<tr><th scope="row">JMS9-P8c</th><td>GSM1545545_JMS9-P8c.txt </td><td>LP                      </td><td>  0                     </td><td>1                       </td><td>L008                    </td></tr>
</tbody>
</table>
</dd>
	<dt>$counts</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">10_6_5_11</th><th scope="col">9_6_5_11</th><th scope="col">purep53</th><th scope="col">JMS8-2</th><th scope="col">JMS8-3</th><th scope="col">JMS8-4</th><th scope="col">JMS8-5</th><th scope="col">JMS9-P7c</th><th scope="col">JMS9-P8c</th></tr></thead>
<tbody>
	<tr><th scope="row">497097</th><td>1  </td><td>2  </td><td>342</td><td>526</td><td>3  </td><td>3  </td><td>535</td><td>2  </td><td>0  </td></tr>
	<tr><th scope="row">100503874</th><td>0  </td><td>0  </td><td>  5</td><td>  6</td><td>0  </td><td>0  </td><td>  5</td><td>0  </td><td>0  </td></tr>
	<tr><th scope="row">100038431</th><td>0  </td><td>0  </td><td>  0</td><td>  0</td><td>0  </td><td>0  </td><td>  1</td><td>0  </td><td>0  </td></tr>
</tbody>
</table>
</dd>
	<dt>$genes</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th scope="col">ENTREZID</th><th scope="col">SYMBOL</th><th scope="col">TXCHROM</th></tr></thead>
<tbody>
	<tr><td>497097   </td><td>Xkr4     </td><td>chr1     </td></tr>
	<tr><td>100503874</td><td>Gm19938  </td><td>NA       </td></tr>
	<tr><td>100038431</td><td>Gm10568  </td><td>NA       </td></tr>
</tbody>
</table>
</dd>
</dl>



Genes that are not expressed at a biologically meaningful level in any condition
should be discarded to reduce the subset of genes to those that are of interest,
and to reduce the number of tests carried out downstream when looking at
differential expression. Upon examination of log-CPM values, it can be seen that
a large proportion of genes within each sample is unexpressed or lowly-expressed
(shown in panel A of the next figure). Using a nominal CPM value of 1 (which is
equivalent to a log-CPM value of 0) genes are deemed to be expressed if their
expression is above this threshold, and unexpressed otherwise. Genes must be
expressed in at least one group (or in at least three samples across the entire
experiment) to be kept for downstream analysis.


Although any sensible value can be used as the expression cutoff, typically a
CPM value of 1 is used in our analyses as it separates expressed genes from
unexpressed genes well for most datasets. Here, a CPM value of 1 means that a
gene is expressed if it has at least 20 counts in the sample with the lowest
sequencing depth (JMS9-P8c, library size approx. 20 million) or at least 76
counts in the sample with the greatest sequencing depth (JMS8-3, library size
approx. 76 million). If sequence reads are summarised by exons rather than genes
and/or experiments have low sequencing depth, a lower CPM cutoff may be
considered.

<br>
<font color ='#00bcd4'> In [22]: </font>

{% highlight R %}
#When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample (cf. the y$samples data.frame)
#will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs, ,keep.lib.sizes=FALSE]
{% endhighlight %}

Using this criterion, the number of genes is reduced to approximately half the
number that we started with (14,165 genes, panel B of the next figure). Note
that subsetting the entire DGEList-object removes both the counts as well as the
associated gene information. Code to produce the figure is given below.

<br>
<font color ='#00bcd4'> In [23]: </font>

{% highlight R %}
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
#lcpm stores the original log data
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
 den <- density(lcpm[,i])
 lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

#compute the log of the filtered data
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
   den <- density(lcpm[,i])
   lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_36_0.png)


The density of log-CPM values for raw pre-filtered data (A) and post-filtered
data (B) are shown for each sample. Dotted vertical lines mark the log-CPM of
zero threshold (equivalent to a CPM value of 1) used in the filtering step.


### Normalising gene expression distributions

During the sample preparation or sequencing process, external factors that are
not of biological interest can affect the expression of individual samples. For
example, samples processed in the first batch of an experiment can have higher
expression overall when compared to samples processed in a second batch. It is
assumed that all samples should have a similar range and distribution of
expression values. Normalisation is required to ensure that the expression
distributions of each sample are similar across the entire experiment.

Any plot showing the per sample expression distributions, such as a density or
boxplot, is useful in determining whether any samples are dissimilar to others.
Distributions of log-CPM values are similar throughout all samples within this
dataset (panel B of the figure above).

Nonetheless, normalisation by the method of trimmed mean of M-values (TMM)
(Robinson and Oshlack 2010) is performed using the `calcNormFactors` function in
edgeR. The normalisation factors calculated here are used as a scaling factor
for the library sizes. When working with DGEList-objects, these normalisation
factors are automatically stored in `x$samples$norm.factors`. For this dataset
the effect of TMM-normalisation is mild, as evident in the magnitude of the
scaling factors, which are all relatively close to 1.

<br>
<font color ='#00bcd4'> In [24]: </font>

{% highlight R %}
x <- calcNormFactors(x, method = "TMM")

x$samples$norm.factors
{% endhighlight %}


<ol class="list-inline">
	<li>0.895730895914575</li>
	<li>1.03491959450621</li>
	<li>1.0439552072696</li>
	<li>1.04050399903738</li>
	<li>1.03235994194938</li>
	<li>0.922342409820323</li>
	<li>0.983660280424014</li>
	<li>1.08273812106367</li>
	<li>0.979260713843751</li>
</ol>



To give a better visual representation of the effects of normalisation, the data
was duplicated then adjusted so that the counts of the first sample are reduced
to 5% of their original values, and in the second sample they are inflated to be
5-times larger.

<br>
<font color ='#00bcd4'> In [25]: </font>

{% highlight R %}
#this is a fake to see the effect of normalization per library size!
#Just see the effect on the first two samples (the other are normalized already)!
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
{% endhighlight %}

The figure below shows the expression distribution of samples for unnormalised
and normalised data, where distributions are noticeably different pre-
normalisation and are similar post-normalisation. Here the first sample has a
small TMM scaling factor of 0.05, whereas the second sample has a large scaling
factor of 6.13 – neither values are close to 1.

<br>
<font color ='#00bcd4'> In [26]: </font>

{% highlight R %}
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
{% endhighlight %}


<ol class="list-inline">
	<li>0.0547222259540704</li>
	<li>6.13059439546747</li>
	<li>1.22927354811535</li>
	<li>1.17051887408206</li>
	<li>1.21487709006519</li>
	<li>1.05622967906675</li>
	<li>1.14587663386915</li>
	<li>1.26129349554649</li>
	<li>1.11702263916712</li>
</ol>




![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_43_1.png)


Example data: Boxplots of log-CPM values showing expression distributions for
unnormalised data (A) and normalised data (B) for each sample in the modified
dataset where the counts in samples 1 and 2 have been scaled to 5% and 500% of
their original values respectively.

### Unsupervised clustering of samples


One of the most important exploratory plots to examine for gene expression
analyses is the multi-dimensional scaling (MDS) plot, or similar. The plot shows
similarities and dissimilarities between samples in an unsupervised manner so
that one can have an idea of the extent to which differential expression can be
detected before carrying out formal tests. Ideally, samples would cluster well
within the primary condition of interest, and any sample straying far from its
group could be identified and followed up for sources of error or extra
variation. If present, technical replicates should lie very close to one
another.

Such a plot can be made in `limma` using the `plotMDS` function. The first
dimension represents the leading-fold-change that best separates samples and
explains the largest proportion of variation in the data, with subsequent
dimensions having a smaller effect and being orthogonal to the ones before it.
When experimental design involves multiple factors, it is recommended that each
factor is examined over several dimensions. If samples cluster by a given factor
in any of these dimensions, it suggests that the factor contributes to
expression differences and is worth including in the linear modelling. On the
other hand, factors that show little or no effect may be left out of downstream
analysis.

In this dataset, samples can be seen to cluster well within experimental groups
over dimension 1 and 2, and then separate by sequencing lane (sample batch) over
dimension 3 (shown in the plot below). Keeping in mind that the first dimension
explains the largest proportion of variation in the data, notice that the range
of values over the dimensions become smaller as we move to higher dimensions.
Whilst all samples cluster by groups, the largest transcriptional difference is
observed between basal and LP, and basal and ML over dimension 1. For this
reason, it is expected that pairwise comparisons between cell populations will
result in a greater number of DE genes for comparisons involving basal samples,
and relatively small numbers of DE genes when comparing ML to LP. In other
datasets, samples that do not cluster by their groups of interest may also show
little or no evidence of differential expression in the downstream analysis.
To create the MDS plots, different colour groupings are assigned to factors of
interest. Dimensions 1 and 2 are examined using the color grouping defined by
cell types.
Dimensions 3 and 4 are examined using the colour grouping defined by sequencing
lanes (batch).

<br>
<font color ='#00bcd4'> In [27]: </font>

{% highlight R %}
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_47_0.png)


MDS plots of log-CPM values over dimensions 1 and 2 with samples coloured and
labeled by sample groups (A) and over dimensions 3 and 4 with samples coloured
and labeled by sequencing lane (B). Distances on the plot correspond to the
leading fold-change, which is the average (root-mean-square) log2-fold-change
for the 500 genes most divergent between each pair of samples by default.

Alternatively, the Glimma package offers the convenience of an interactive MDS
plot where multiple dimensions can be explored. The glMDSPlot function generates
an html page (that is opened in a browser if launch=TRUE) with an MDS plot in
the left panel and a barplot showing the proportion of variation explained by
each dimension in the right panel. Clicking on the bars of the bar plot changes
the pair of dimensions plotted in the MDS plot, and hovering over the individual
points reveals the sample label. The colour scheme can be changed as well to
highlight cell population or sequencing lane (batch). An interactive MDS plot of
this dataset can be found at
<http://bioinf.wehi.edu.au/folders/limmaWorkflow/glimma-plots/MDS-Plot.html>.

<br>
<font color ='#00bcd4'> In [28]: </font>

{% highlight R %}
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"),
          groups=x$samples[,c(2,5)], launch=FALSE)
{% endhighlight %}

## Differential expression analysis
### Creating a design

In this study, it is of interest to see which genes are expressed at different
levels between the three cell populations profiled. In our analysis, linear
models are fitted to the data with the assumption that the underlying data is
normally distributed. To get started, a design matrix is set up with both the
cell population and sequencing lane (batch) information.

<br>
<font color ='#00bcd4'> In [29]: </font>

{% highlight R %}
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))

design
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">Basal</th><th scope="col">LP</th><th scope="col">ML</th><th scope="col">laneL006</th><th scope="col">laneL008</th></tr></thead>
<tbody>
	<tr><th scope="row">1</th><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope="row">2</th><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope="row">3</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope="row">4</th><td>1</td><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope="row">5</th><td>0</td><td>0</td><td>1</td><td>1</td><td>0</td></tr>
	<tr><th scope="row">6</th><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope="row">7</th><td>1</td><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope="row">8</th><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td></tr>
	<tr><th scope="row">9</th><td>0</td><td>1</td><td>0</td><td>0</td><td>1</td></tr>
</tbody>
</table>



For a given experiment, there are usually several equivalent ways to set up an
appropriate design matrix. For example, `~0+group+lane` removes the intercept
from the first factor, `group`, but an intercept remains in the second factor
`lane`. Alternatively, `~group+lane` could be used to keep the intercepts in
both group and lane. Understanding how to interpret the coefficients estimated
in a given model is key here. We choose the first model for our analysis, as
setting up model contrasts is more straight forward in the absence of an
intercept for `group`. Contrasts for pairwise comparisons between cell
populations are set up in limma using the `makeContrasts` function.

<br>
<font color ='#00bcd4'> In [30]: </font>

{% highlight R %}
contr.matrix <- makeContrasts(
   BasalvsLP = Basal-LP,
   BasalvsML = Basal - ML,
   LPvsML = LP - ML,
   levels = colnames(design))
contr.matrix
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">BasalvsLP</th><th scope="col">BasalvsML</th><th scope="col">LPvsML</th></tr></thead>
<tbody>
	<tr><th scope="row">Basal</th><td> 1</td><td> 1</td><td> 0</td></tr>
	<tr><th scope="row">LP</th><td>-1</td><td> 0</td><td> 1</td></tr>
	<tr><th scope="row">ML</th><td> 0</td><td>-1</td><td>-1</td></tr>
	<tr><th scope="row">laneL006</th><td> 0</td><td> 0</td><td> 0</td></tr>
	<tr><th scope="row">laneL008</th><td> 0</td><td> 0</td><td> 0</td></tr>
</tbody>
</table>



A key strength of limma’s linear modelling approach, is the ability accommodate
arbitrary experimental complexity. Simple designs, such as the one in this
workflow, with cell type and batch, through to more complicated factorial
designs and models with interaction terms can be handled relatively easily.
Where experimental or technical effects can be modelled using a random effect,
another possibility in limma is to estimate correlations using
`duplicateCorrelation` by specifying a `block` argument for both this function
and in the `lmFit` linear modelling step.

### Removing heteroscedascity from count data

It has been shown that for RNA-seq count data, the variance is not independent
of the mean (Law et al. 2014) – this is true of raw counts or when transformed
to log-CPM values. Methods that model counts using a Negative Binomial
distribution assume a quadratic mean-variance relationship. In limma, linear
modelling is carried out on the log-CPM values which are assumed to be normally
distributed and the mean-variance relationship is accommodated using precision
weights calculated by the `voom` function.

When operating on a DGEList-object, `voom` converts raw counts to log-CPM values
by automatically extracting library sizes and normalisation factors from x
itself. Additional normalisation to log-CPM values can be specified within voom
using the `normalize.method` argument.

The mean-variance relationship of log-CPM values for this dataset is shown in
the left-hand panel of the next figure. Typically, the voom-plot shows a
decreasing trend between the means and variances resulting from a combination of
technical variation in the sequencing experiment and biological variation
amongst the replicate samples from different cell populations. Experiments with
high biological variation usually result in flatter trends, where variance
values plateau at high expression values. Experiments with low biological
variation tend to result in sharp decreasing trends.

Moreover, the voom-plot provides a visual check on the level of filtering
performed upstream. If filtering of lowly-expressed genes is insufficient, a
drop in variance levels can be observed at the low end of the expression scale
due to very small counts. If this is observed, one should return to the earlier
filtering step and increase the expression threshold applied to the dataset.

Where sample-level variation is evident from earlier inspections of the MDS
plot, the `voomWithQualityWeights` function can be used to simultaneously
incorporate sample-level weights together with the abundance dependent weights
estimated by `voom` (Liu et al. 2015). For an example of this, see Liu et al.
(2016) (Liu et al. 2016).

<br>
<font color ='#00bcd4'> In [31]: </font>

{% highlight R %}
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
#v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_56_1.png)


Note that the other data frames stored within the `DGEList-object` that contain
gene- and sample-level information, are retained in the EList-object `v` created
by `voom`. The `v$genes` data frame is equivalent to `x$genes`, `v$targets` is
equivalent to `x$samples`, and the expression values stored in v$E is analogous
to `x$counts`, albeit on a transformed scale. In addition to this, the `voom`
EList-object has a matrix of precision weights `v$weights` and stores the design
matrix in `v$design`.

### Fitting linear models for comparisons of interest

Linear modelling in limma is carried out using the `lmFit` and `contrasts.fit`
functions originally written for application to microarrays. The functions can
be used for both microarray and RNA-seq data and fit a separate model to the
expression values for each gene. Next, empirical Bayes moderation is carried out
by borrowing information across all the genes to obtain more precise estimates
of gene-wise variability (Smyth 2004). The model’s residual variances are
plotted against average expression values in the next figure. It can be seen
from this plot that the variance is no longer dependent on the mean expression
level.

### Examining the number of DE genes

For a quick look at differential expression levels, the number of significantly
up- and down-regulated genes can be summarised in a table. Significance is
defined using an adjusted p-value cutoff that is set at 5% by default. For the
comparison between expression levels in basal and LP, 4,127 genes are found to
be down-regulated in basal relative to LP and 4,298 genes are up-regulated in
basal relative to LP – a total of 8,425 DE genes. A total of 8,510 DE genes are
found between basal and ML (4,338 down- and 4,172 up-regulated genes), and a
total of 5,340 DE genes are found between LP and ML (2,895 down- and 2,445 up-
regulated). The larger numbers of DE genes observed for comparisons involving
the basal population are consistent with our observations from the MDS plots.

<br>
<font color ='#00bcd4'> In [32]: </font>

{% highlight R %}
summary(decideTests(efit))
{% endhighlight %}


       BasalvsLP BasalvsML LPvsML
    -1      4127      4338   2895
    0       5740      5655   8825
    1       4298      4172   2445


Some studies require more than an adjusted p-value cut-off. For a stricter
definition on significance, one may require log-fold-changes (log-FCs) to be
above a minimum value. The treat method (McCarthy and Smyth 2009) can be used to
calculate p-values from empirical Bayes moderated t-statistics with a minimum
log-FC requirement. The number of differentially expressed genes are reduced to
a total of 3,135 DE genes for basal versus LP, 3,270 DE genes for basal versus
ML, and 385 DE genes for LP versus ML when testing requires genes to have a log-
FC that is significantly greater than 1 (equivalent to a 2-fold difference
between cell types on the original scale).

<br>
<font color ='#00bcd4'> In [33]: </font>

{% highlight R %}
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
{% endhighlight %}


       BasalvsLP BasalvsML LPvsML
    -1      1417      1512    203
    0      11030     10895  13780
    1       1718      1758    182


Genes that are DE in multiple comparisons can be extracted using the results
from `decideTests`, where 0s represent genes that are not DE, 1s represent genes
that are up-regulated, and -1s represent genes that are down-regulated. A total
of 2,409 genes are DE in both basal versus LP and basal versus ML, twenty of
which are listed below. The `write.fit` function can be used to extract and
write results for all three comparisons to a single output file.

<br>
<font color ='#00bcd4'> In [34]: </font>

{% highlight R %}
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
{% endhighlight %}


2409


<br>
<font color ='#00bcd4'> In [35]: </font>

{% highlight R %}
head(tfit$genes$SYMBOL[de.common], n=20)
{% endhighlight %}


<ol class="list-inline">
	<li>'Xkr4'</li>
	<li>'Rgs20'</li>
	<li>'Cpa6'</li>
	<li>'Sulf1'</li>
	<li>'Eya1'</li>
	<li>'Msc'</li>
	<li>'Sbspon'</li>
	<li>'Pi15'</li>
	<li>'Crispld1'</li>
	<li>'Kcnq5'</li>
	<li>'Ptpn18'</li>
	<li>'Arhgef4'</li>
	<li>'2010300C02Rik'</li>
	<li>'Aff3'</li>
	<li>'Npas2'</li>
	<li>'Tbc1d8'</li>
	<li>'Creg2'</li>
	<li>'Il1r1'</li>
	<li>'Il18r1'</li>
	<li>'Il18rap'</li>
</ol>



<br>
<font color ='#00bcd4'> In [36]: </font>

{% highlight R %}
vennDiagram(dt[,1:2], circle.col=c("mediumpurple", "midnightblue"))
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_64_0.png)


<br>
<font color ='#00bcd4'> In [37]: </font>

{% highlight R %}
vennDiagram(dt[,1:3], circle.col=c("mediumpurple", "midnightblue"))
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_65_0.png)


Venn diagram showing the number of genes DE in the comparison between basal
versus LP only (left), basal versus ML only (right), and the number of genes
that are DE in both comparisons (center). The number of genes that are not DE in
either comparison are marked in the bottom-right.

<br>
<font color ='#00bcd4'> In [38]: </font>

{% highlight R %}
write.fit(tfit, dt, file="results.txt")
{% endhighlight %}

### Examining individual DE genes from top to bottom

The top DE genes can be listed using `topTreat` for results using treat (or
`topTable` for results using `eBayes`). By default `topTreat` arranges genes
from smallest to largest adjusted p-value with associated gene information, log-
FC, average log-CPM, moderated t-statistic, raw and adjusted p-value for each
gene. The number of top genes displayed can be specified, where `n=Inf` includes
all genes. Genes Cldn7 and Rasef are amongst the top DE genes for both basal
versus LP and basal versus ML.

<br>
<font color ='#00bcd4'> In [39]: </font>

{% highlight R %}
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">ENTREZID</th><th scope="col">SYMBOL</th><th scope="col">TXCHROM</th><th scope="col">logFC</th><th scope="col">AveExpr</th><th scope="col">t</th><th scope="col">P.Value</th><th scope="col">adj.P.Val</th></tr></thead>
<tbody>
	<tr><th scope="row">12759</th><td>12759       </td><td>Clu         </td><td>chr14       </td><td>-5.442877   </td><td>8.857907    </td><td>-33.44429   </td><td>3.990899e-10</td><td>2.703871e-06</td></tr>
	<tr><th scope="row">53624</th><td>53624       </td><td>Cldn7       </td><td>chr11       </td><td>-5.514605   </td><td>6.296762    </td><td>-32.94533   </td><td>4.503694e-10</td><td>2.703871e-06</td></tr>
	<tr><th scope="row">242505</th><td>242505      </td><td>Rasef       </td><td>chr4        </td><td>-5.921741   </td><td>5.119585    </td><td>-31.77625   </td><td>6.063249e-10</td><td>2.703871e-06</td></tr>
	<tr><th scope="row">67451</th><td>67451       </td><td>Pkp2        </td><td>chr16       </td><td>-5.724823   </td><td>4.420495    </td><td>-30.65370   </td><td>8.010456e-10</td><td>2.703871e-06</td></tr>
	<tr><th scope="row">228543</th><td>228543      </td><td>Rhov        </td><td>chr2        </td><td>-6.253427   </td><td>5.486640    </td><td>-29.46244   </td><td>1.112729e-09</td><td>2.703871e-06</td></tr>
	<tr><th scope="row">70350</th><td>70350       </td><td>Basp1       </td><td>chr15       </td><td>-6.073297   </td><td>5.248349    </td><td>-28.64890   </td><td>1.380545e-09</td><td>2.703871e-06</td></tr>
</tbody>
</table>



### Useful graphical representations of differential expression results

To summarise results for all genes visually, mean-difference plots, which
display log-FCs from the linear model fit against the average log-CPM values can
be generated using the `plotMD` function, with the differentially expressed
genes highlighted.

<br>
<font color ='#00bcd4'> In [40]: </font>

{% highlight R %}
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1],
       xlim=c(-8,13), col = c('mediumpurple', 'orange'))
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_71_0.png)


*Glimma* extends this functionality by providing an interactive mean-difference
plot via the `glMDPlot` function. The output of this function is an html page,
with summarised results in the left panel (similar to what is output by
`plotMD`), and the `log-CPM` values from individual samples in the right panel,
with a table of results below the plots. This interactive display allows the
user to search for particular genes based on their Gene symbol, which is not
possible in a static *R* plot. The `glMDPlot` function is not limited to mean-
difference plots, with a default version allowing a data frame to be passed with
the user able to select the columns of interest to plot in the left panel
(interactive plot
[here](http://www.bioconductor.org/help/workflows/RNAseq123/glimma-plots/MD-
Plot.html)).

<br>
<font color ='#00bcd4'> In [41]: </font>

{% highlight R %}
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         id.column="ENTREZID", counts=x$counts, groups=group, launch=FALSE)
{% endhighlight %}

Interactive mean-difference plot generated using Glimma. Summary data (log-FCs
versus log-CPM values) are shown in the left panel which is linked to the
individual values per sample for a selected gene in the right panel. A table of
results is also displayed below these figures, along with a search bar to allow
users to look up a particular gene using the annotation information available,
e.g. the Gene symbol identifier Clu.

The mean-difference plot generated by the command above is available online (see
<http://bioinf.wehi.edu.au/folders/limmaWorkflow/glimma-plots/MD-Plot.html>).
The interactivity provided by the Glimma package allows additional information
to be presented in a single graphical window. Glimma is implemented in R and
Javascript, with the R code generating the data which is converted into graphics
using the Javascript library D3 (<https://d3js.org>), with the Bootstrap library
handling layouts and Datatables generating the interactive searchable tables.
This allows plots to be viewed in any modern browser, which is convenient for
including them as linked files from an Rmarkdown report of the analysis.

Plots shown previously include either all of the genes that are expressed in any
one condition (such as the Venn diagram of common DE genes or mean-difference
plot) or look at genes individually (log-CPM values shown in right panel of the
interactive mean-difference plot). Heatmaps allow users to look at the
expression of a subset of genes. This can be give useful insight into the
expression of individual groups and samples without losing perspective of the
overall study when focusing on individual genes, or losing resolution when
examining patterns averaged over thousands of genes at the same time.

A heatmap is created for the top 100 DE genes (as ranked by adjusted `p-value`)
from the basal versus LP contrast using the `heatmap.2` function from the
_gplots_ package. The heatmap correctly clusters samples into cell type and
rearranges the order of genes to form blocks of similar expression. From the
heatmap, we observe that the expression of ML and LP samples are very similar
for the top 100 DE genes between basal and LP.

<br>
<font color ='#00bcd4'> In [42]: </font>

{% highlight R %}
library(gplots)

basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"orange","white","darkslateblue")
heatmap.2(v$E[i,], scale="row",
   labRow=v$genes$SYMBOL[i], labCol=group,
   col=mycol, trace="none", density.info="none",
   margin=c(8,6), lhei=c(2,10), dendrogram="column")
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_75_1.png)


Heatmap of log-CPM values for top 100 genes DE in basal versus LP. Expression
across each gene (or row) have been scaled so that mean expression is zero and
standard deviation is one. Samples with relatively high expression of a given
gene are marked in red and samples with relatively low expression are marked in
blue. Lighter shades and white represent genes with intermediate expression
levels. Samples and genes have been reordered by the method of hierarchical
clustering. A dendrogram is shown for the sample clustering.

## Gene testing with camera

We finish off this analysis with some gene set testing by applying the camera
method (Wu and Smyth 2012) to the c2 gene signatures from Broad Institute’s
MSigDB c2 collection (Subramanian et al. 2005) that have been adapted for mouse
and are available as Rdata objects from
<http://bioinf.wehi.edu.au/software/MSigDB/>. Other useful gene sets derived
from MSigDB for both human and mouse, such as the hallmark gene sets, are also
available from this site. C2 gene sets have been curated from online databases,
publications and domain experts, and hallmark gene sets are selected from MSigDB
to have well-defined biological states or processes.

<br>
<font color ='#00bcd4'> In [43]: </font>

{% highlight R %}
load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
idx <- ids2indices(Mm.c2,id=rownames(v))
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.BasalvsLP,5)
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">NGenes</th><th scope="col">Direction</th><th scope="col">PValue</th><th scope="col">FDR</th></tr></thead>
<tbody>
	<tr><th scope="row">LIM_MAMMARY_STEM_CELL_UP</th><td>739         </td><td>Up          </td><td>1.134757e-18</td><td>5.360590e-15</td></tr>
	<tr><th scope="row">LIM_MAMMARY_STEM_CELL_DN</th><td>630         </td><td>Down        </td><td>1.569957e-15</td><td>3.708238e-12</td></tr>
	<tr><th scope="row">ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER</th><td>163         </td><td>Up          </td><td>1.437987e-13</td><td>2.264351e-10</td></tr>
	<tr><th scope="row">SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP</th><td>183         </td><td>Up          </td><td>2.181862e-13</td><td>2.576779e-10</td></tr>
	<tr><th scope="row">LIM_MAMMARY_LUMINAL_PROGENITOR_UP</th><td> 87         </td><td>Down        </td><td>6.734613e-13</td><td>6.362863e-10</td></tr>
</tbody>
</table>



<br>
<font color ='#00bcd4'> In [44]: </font>

{% highlight R %}
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.BasalvsML,5)
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">NGenes</th><th scope="col">Direction</th><th scope="col">PValue</th><th scope="col">FDR</th></tr></thead>
<tbody>
	<tr><th scope="row">LIM_MAMMARY_STEM_CELL_UP</th><td>739         </td><td>Up          </td><td>5.090937e-23</td><td>2.404959e-19</td></tr>
	<tr><th scope="row">LIM_MAMMARY_STEM_CELL_DN</th><td>630         </td><td>Down        </td><td>5.132446e-19</td><td>1.212284e-15</td></tr>
	<tr><th scope="row">LIM_MAMMARY_LUMINAL_MATURE_DN</th><td>166         </td><td>Up          </td><td>8.875174e-16</td><td>1.397544e-12</td></tr>
	<tr><th scope="row">LIM_MAMMARY_LUMINAL_MATURE_UP</th><td>180         </td><td>Down        </td><td>6.287301e-13</td><td>7.425303e-10</td></tr>
	<tr><th scope="row">ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER</th><td>163         </td><td>Up          </td><td>1.684323e-12</td><td>1.591348e-09</td></tr>
</tbody>
</table>



<br>
<font color ='#00bcd4'> In [45]: </font>

{% highlight R %}
cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.LPvsML,5)
{% endhighlight %}


<table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">NGenes</th><th scope="col">Direction</th><th scope="col">PValue</th><th scope="col">FDR</th></tr></thead>
<tbody>
	<tr><th scope="row">LIM_MAMMARY_LUMINAL_MATURE_UP</th><td>180         </td><td>Down        </td><td>8.497295e-14</td><td>3.401020e-10</td></tr>
	<tr><th scope="row">LIM_MAMMARY_LUMINAL_MATURE_DN</th><td>166         </td><td>Up          </td><td>1.439890e-13</td><td>3.401020e-10</td></tr>
	<tr><th scope="row">LIM_MAMMARY_LUMINAL_PROGENITOR_UP</th><td> 87         </td><td>Up          </td><td>3.840915e-11</td><td>6.048160e-08</td></tr>
	<tr><th scope="row">REACTOME_RESPIRATORY_ELECTRON_TRANSPORT</th><td> 91         </td><td>Down        </td><td>2.655349e-08</td><td>3.135967e-05</td></tr>
	<tr><th scope="row">NABA_CORE_MATRISOME</th><td>222         </td><td>Up          </td><td>4.430361e-08</td><td>4.185805e-05</td></tr>
</tbody>
</table>



The `camera` function performs a competitive test to assess whether the genes in
a given set are highly ranked in terms of differential expression relative to
genes that are not in the set. It uses limma’s linear model framework, taking
both the design matrix and contrast matrix (if present) and accommodates the
observational-level weights from voom in the testing procedure. After adjusting
the variance of the resulting gene set test statistic by a variance inflation
factor, that depends on the gene-wise correlation (which is set to 0.01 by
default, but can be estimated from the data) and the size of the set, a p-value
is returned and adjusted for multiple testing.


This experiment is the RNA-seq equivalent of Lim et al. (2010) (Lim et al.
2010), who used Illumina microarrays to profile the same sorted cell
populations, so it is reassuring to see the gene signatures from this earlier
publication coming up at the top of the list for each contrast. We make a
barcodeplot of the Lim et al. (2010) Mature Luminal gene sets (Up and Down) in
the LP versus ML contrast. Note that these sets go in the opposite direction in
our dataset due to our parameterization which compares LP against ML rather than
the other way around (if the contrast were reversed, the directions would be
consistent).

<br>
<font color ='#00bcd4'> In [46]: </font>

{% highlight R %}
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP,
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/rna-sequencing_files/rna-sequencing_82_0.png)


Barcode plot of `LIM_MAMMARY_LUMINAL_MATURE_UP` (red bars, top of plot) and
`LIM_MAMMARY_LUMINAL_MATURE_DN` (blue bars, bottom of plot) gene sets in the LP
versus ML contrast with an enrichment line for each set that shows the relative
enrichment of the vertical bars in each part of the plot. The experiment of Lim
et al. (Lim et al. 2010) is very similar to the current one, with the same
sorting strategy used to obtain the different cell populations, except that
microarrays were used instead of RNA-seq to profile gene expression. Note that
the inverse correlation (the up gene set is down and the down gene set is up) is
a result of the way the contrast has been set up (LP versus ML) – if reversed,
the directionality would agree.


Other gene set tests are available in limma, such as the self-contained tests by
mroast (Wu et al. 2010). Whilst camera is ideal for testing a large database of
gene sets and observing which of them rank highly relative to others (as shown
above), self-contained tests are better for focused testing of one or a few
specifically chosen sets to see if they are DE in their own right. In other
words, camera is more appropriate when ‘fishing’ for gene sets of interest,
whereas mroast tests sets that are already of interest for significance.

---

<a href="{{site.url}}{{site.baseurl}}{% link Day3/MDPlot.md %}" class="button button-default" >
<i class="fa fa-bar-chart" aria-hidden="true"></i> Interactive MD plot
</a>


<a href="{{site.url}}{{site.baseurl}}{% link Day3/MDSPlot.md %}" class="button" >
<i class="fa fa-bar-chart" aria-hidden="true"></i>Interactive MDS plot
</a>



---
