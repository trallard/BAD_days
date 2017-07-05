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

## Loading the required packages

<br>
<font color ='#00bcd4'> In [1]: </font>

{% highlight R %}
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus);
{% endhighlight %}

## Getting the data

The dataset `GSE63310_RAW.tar` will be directly downloaded from
<https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file>. We will
then extract the relevant files.

<br>
<font color ='#00bcd4'> In [2]: </font>

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
<font color ='#00bcd4'> In [3]: </font>

{% highlight R %}
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",
   "GSM1545545_JMS9-P8c.txt")

read.delim(files[1], nrow=5)
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



Instead of readig each file separately and combine posteriorly the package
`edgeR` allows us to do this in one single step using the function `readDGE`.

The resulting object will be a matrix of counts with 27,179 rows associated with
unique Entrez gene identifiers (IDs) and nine columns associated with the
individual samples in the experiment.

<br>
<font color ='#00bcd4'> In [4]: </font>

{% highlight R %}
x <- readDGE(files, columns=c(1,3))
class(x)
{% endhighlight %}


'DGEList'


<br>
<font color ='#00bcd4'> In [5]: </font>

{% highlight R %}
# checking the dimension of x
dim(x)
{% endhighlight %}


<ol class="list-inline">
	<li>27179</li>
	<li>9</li>
</ol>



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
<font color ='#00bcd4'> In [6]: </font>

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
<font color ='#00bcd4'> In [7]: </font>

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
<font color ='#00bcd4'> In [8]: </font>

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
<font color ='#00bcd4'> In [9]: </font>

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
<font color ='#00bcd4'> In [10]: </font>

{% highlight R %}
x$genes <- genes
x
{% endhighlight %}


<dl>
	<dt>$samples</dt>
		<dd><table class="table-responsive table-striped">
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
</dd>
	<dt>$counts</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">10_6_5_11</th><th scope="col">9_6_5_11</th><th scope="col">purep53</th><th scope="col">JMS8-2</th><th scope="col">JMS8-3</th><th scope="col">JMS8-4</th><th scope="col">JMS8-5</th><th scope="col">JMS9-P7c</th><th scope="col">JMS9-P8c</th></tr></thead>
<tbody>
	<tr><th scope="row">497097</th><td>   1</td><td>   2</td><td> 342</td><td> 526</td><td>   3</td><td>   3</td><td> 535</td><td>   2</td><td>   0</td></tr>
	<tr><th scope="row">100503874</th><td>   0</td><td>   0</td><td>   5</td><td>   6</td><td>   0</td><td>   0</td><td>   5</td><td>   0</td><td>   0</td></tr>
	<tr><th scope="row">100038431</th><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   1</td><td>   0</td><td>   0</td></tr>
	<tr><th scope="row">19888</th><td>   0</td><td>   1</td><td>   0</td><td>   0</td><td>  17</td><td>   2</td><td>   0</td><td>   1</td><td>   0</td></tr>
	<tr><th scope="row">20671</th><td>   1</td><td>   1</td><td>  76</td><td>  40</td><td>  33</td><td>  14</td><td>  98</td><td>  18</td><td>   8</td></tr>
	<tr><th scope="row">27395</th><td> 431</td><td> 771</td><td>1368</td><td>1268</td><td>1564</td><td> 769</td><td> 818</td><td> 468</td><td> 342</td></tr>
	<tr><th scope="row">18777</th><td> 768</td><td>1722</td><td>2517</td><td>1923</td><td>3865</td><td>1888</td><td>1830</td><td>1246</td><td> 693</td></tr>
	<tr><th scope="row">100503730</th><td>   4</td><td>   8</td><td>   6</td><td>   2</td><td>  11</td><td>  11</td><td>   3</td><td>   9</td><td>   2</td></tr>
	<tr><th scope="row">21399</th><td> 810</td><td> 977</td><td>2472</td><td>1870</td><td>2251</td><td>1716</td><td>1932</td><td> 756</td><td> 619</td></tr>
	<tr><th scope="row">58175</th><td> 452</td><td> 358</td><td>  17</td><td>  14</td><td> 622</td><td> 571</td><td>  12</td><td> 203</td><td> 224</td></tr>
	<tr><th scope="row">108664</th><td>1716</td><td>2678</td><td>2097</td><td>2071</td><td>5499</td><td>3630</td><td>1731</td><td>1715</td><td>1251</td></tr>
	<tr><th scope="row">18387</th><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope="row">226304</th><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope="row">12421</th><td>3451</td><td>2699</td><td>3399</td><td>2716</td><td>5233</td><td>6280</td><td>3647</td><td>1866</td><td>2122</td></tr>
	<tr><th scope="row">620393</th><td>   0</td><td>   0</td><td>   2</td><td>   3</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope="row">240690</th><td>   0</td><td>   0</td><td>   0</td><td>   1</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope="row">319263</th><td>2026</td><td>2033</td><td>3920</td><td>2715</td><td>3873</td><td>3688</td><td>3593</td><td>1609</td><td>1348</td></tr>
	<tr><th scope="row">71096</th><td>   0</td><td>   0</td><td>   1</td><td>   5</td><td>   0</td><td>   0</td><td>   4</td><td>   0</td><td>   0</td></tr>
	<tr><th scope="row">59014</th><td> 956</td><td> 985</td><td>5497</td><td>4214</td><td>3462</td><td>2933</td><td>5336</td><td> 649</td><td> 731</td></tr>
	<tr><th scope="row">76187</th><td>  54</td><td>  76</td><td>  46</td><td>  32</td><td> 148</td><td> 126</td><td>  59</td><td>  44</td><td>  27</td></tr>
	<tr><th scope="row">72481</th><td>  16</td><td>  28</td><td>  12</td><td>  12</td><td>  67</td><td>  24</td><td>  14</td><td>  18</td><td>  11</td></tr>
	<tr><th scope="row">76982</th><td>  21</td><td>   9</td><td>  21</td><td>  44</td><td>  35</td><td>  34</td><td>  29</td><td>   8</td><td>  15</td></tr>
	<tr><th scope="row">17864</th><td> 194</td><td> 166</td><td>1007</td><td> 862</td><td> 632</td><td> 594</td><td> 901</td><td> 100</td><td> 131</td></tr>
	<tr><th scope="row">70675</th><td>2390</td><td>2012</td><td>3500</td><td>3020</td><td>5412</td><td>4404</td><td>4173</td><td>1585</td><td>1522</td></tr>
	<tr><th scope="row">73331</th><td>   3</td><td>   2</td><td>   4</td><td>   1</td><td>   1</td><td>   5</td><td>   2</td><td>   3</td><td>   0</td></tr>
	<tr><th scope="row">170755</th><td>  87</td><td> 323</td><td> 734</td><td> 565</td><td> 721</td><td> 249</td><td> 666</td><td> 267</td><td> 113</td></tr>
	<tr><th scope="row">620986</th><td>  71</td><td> 194</td><td> 383</td><td> 206</td><td> 307</td><td> 110</td><td> 206</td><td> 112</td><td>  50</td></tr>
	<tr><th scope="row">240697</th><td>   0</td><td>   0</td><td>   9</td><td>   8</td><td>   3</td><td>   1</td><td>   7</td><td>   0</td><td>   1</td></tr>
	<tr><th scope="row">73824</th><td> 289</td><td> 487</td><td> 824</td><td> 565</td><td> 709</td><td> 517</td><td> 662</td><td> 335</td><td> 271</td></tr>
	<tr><th scope="row">266793</th><td>   4</td><td>   6</td><td>  54</td><td>  10</td><td>  18</td><td>  19</td><td>  43</td><td>   9</td><td>   6</td></tr>
	<tr><th scope="row">⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope="row">100041168</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100041207</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100861637</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100861988</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100862006</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100862025</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100040991</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100040911</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100862053</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100039905</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100041141</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100041117</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100862083</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100862092</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100862100</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100042201</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100039904</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100861873</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100861881</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100504642</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100041631</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100504702</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100040357</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100861808</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100504460</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100861837</th><td>320 </td><td>472 </td><td> 582</td><td>335 </td><td> 811</td><td> 432</td><td> 500</td><td>200 </td><td> 84 </td></tr>
	<tr><th scope="row">100861924</th><td> 15 </td><td> 27 </td><td> 126</td><td>  7 </td><td>  17</td><td>  22</td><td>  35</td><td>  3 </td><td>  2 </td></tr>
	<tr><th scope="row">170942</th><td>824 </td><td>769 </td><td>1337</td><td>781 </td><td>1266</td><td>1089</td><td>1274</td><td>469 </td><td>379 </td></tr>
	<tr><th scope="row">100861691</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
	<tr><th scope="row">100504472</th><td>  0 </td><td>  0 </td><td>   0</td><td>  0 </td><td>   0</td><td>   0</td><td>   0</td><td>  0 </td><td>  0 </td></tr>
</tbody>
</table>
</dd>
	<dt>$genes</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">ENTREZID</th><th scope="col">SYMBOL</th><th scope="col">TXCHROM</th></tr></thead>
<tbody>
	<tr><th scope="row">1</th><td>497097       </td><td>Xkr4         </td><td>chr1         </td></tr>
	<tr><th scope="row">2</th><td>100503874    </td><td>Gm19938      </td><td>NA           </td></tr>
	<tr><th scope="row">3</th><td>100038431    </td><td>Gm10568      </td><td>NA           </td></tr>
	<tr><th scope="row">4</th><td>19888        </td><td>Rp1          </td><td>chr1         </td></tr>
	<tr><th scope="row">5</th><td>20671        </td><td>Sox17        </td><td>chr1         </td></tr>
	<tr><th scope="row">6</th><td>27395        </td><td>Mrpl15       </td><td>chr1         </td></tr>
	<tr><th scope="row">7</th><td>18777        </td><td>Lypla1       </td><td>chr1         </td></tr>
	<tr><th scope="row">8</th><td>100503730    </td><td>Gm19860      </td><td>NA           </td></tr>
	<tr><th scope="row">9</th><td>21399        </td><td>Tcea1        </td><td>chr1         </td></tr>
	<tr><th scope="row">10</th><td>58175        </td><td>Rgs20        </td><td>chr1         </td></tr>
	<tr><th scope="row">11</th><td>108664       </td><td>Atp6v1h      </td><td>chr1         </td></tr>
	<tr><th scope="row">12</th><td>18387        </td><td>Oprk1        </td><td>chr1         </td></tr>
	<tr><th scope="row">13</th><td>226304       </td><td>Npbwr1       </td><td>chr1         </td></tr>
	<tr><th scope="row">14</th><td>12421        </td><td>Rb1cc1       </td><td>chr1         </td></tr>
	<tr><th scope="row">15</th><td>620393       </td><td>Fam150a      </td><td>chr1         </td></tr>
	<tr><th scope="row">16</th><td>240690       </td><td>St18         </td><td>chr1         </td></tr>
	<tr><th scope="row">17</th><td>319263       </td><td>Pcmtd1       </td><td>chr1         </td></tr>
	<tr><th scope="row">18</th><td>71096        </td><td>Sntg1        </td><td>chr1         </td></tr>
	<tr><th scope="row">19</th><td>59014        </td><td>Rrs1         </td><td>chr1         </td></tr>
	<tr><th scope="row">20</th><td>76187        </td><td>Adhfe1       </td><td>chr1         </td></tr>
	<tr><th scope="row">21</th><td>72481        </td><td>2610203C22Rik</td><td>chr1         </td></tr>
	<tr><th scope="row">22</th><td>76982        </td><td>3110035E14Rik</td><td>chr1         </td></tr>
	<tr><th scope="row">23</th><td>17864        </td><td>Mybl1        </td><td>chr1         </td></tr>
	<tr><th scope="row">24</th><td>70675        </td><td>Vcpip1       </td><td>chr1         </td></tr>
	<tr><th scope="row">25</th><td>73331        </td><td>1700034P13Rik</td><td>chr1         </td></tr>
	<tr><th scope="row">26</th><td>170755       </td><td>Sgk3         </td><td>chr1         </td></tr>
	<tr><th scope="row">27</th><td>620986       </td><td>Gm6195       </td><td>NA           </td></tr>
	<tr><th scope="row">28</th><td>240697       </td><td>Mcmdc2       </td><td>chr1         </td></tr>
	<tr><th scope="row">29</th><td>73824        </td><td>Snhg6        </td><td>chr1         </td></tr>
	<tr><th scope="row">30</th><td>266793       </td><td>Snord87      </td><td>chr1         </td></tr>
	<tr><th scope="row">⋮</th><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope="row">27191</th><td>100041168   </td><td>Gm20863     </td><td>NA          </td></tr>
	<tr><th scope="row">27192</th><td>100041207   </td><td>NA          </td><td>NA          </td></tr>
	<tr><th scope="row">27193</th><td>100861637   </td><td>Gm21095     </td><td>NA          </td></tr>
	<tr><th scope="row">27194</th><td>100861988   </td><td>Gm21380     </td><td>NA          </td></tr>
	<tr><th scope="row">27195</th><td>100862006   </td><td>NA          </td><td>NA          </td></tr>
	<tr><th scope="row">27196</th><td>100862025   </td><td>Gm21409     </td><td>NA          </td></tr>
	<tr><th scope="row">27197</th><td>100040991   </td><td>Gm20856     </td><td>NA          </td></tr>
	<tr><th scope="row">27198</th><td>100040911   </td><td>Gm20854     </td><td>chrY        </td></tr>
	<tr><th scope="row">27199</th><td>100862053   </td><td>Gm21435     </td><td>NA          </td></tr>
	<tr><th scope="row">27200</th><td>100039905   </td><td>Gm20820     </td><td>NA          </td></tr>
	<tr><th scope="row">27201</th><td>100041141   </td><td>Gm20861     </td><td>NA          </td></tr>
	<tr><th scope="row">27202</th><td>100041117   </td><td>Gm20860     </td><td>NA          </td></tr>
	<tr><th scope="row">27203</th><td>100862083   </td><td>Gm21462     </td><td>NA          </td></tr>
	<tr><th scope="row">27204</th><td>100862092   </td><td>Gm21470     </td><td>NA          </td></tr>
	<tr><th scope="row">27205</th><td>100862100   </td><td>Gm21477     </td><td>NA          </td></tr>
	<tr><th scope="row">27206</th><td>100042201   </td><td>Gm20906     </td><td>NA          </td></tr>
	<tr><th scope="row">27207</th><td>100039904   </td><td>Gm20819     </td><td>NA          </td></tr>
	<tr><th scope="row">27208</th><td>100861873   </td><td>Gm21287     </td><td>NA          </td></tr>
	<tr><th scope="row">27209</th><td>100861881   </td><td>Gm21294     </td><td>NA          </td></tr>
	<tr><th scope="row">27210</th><td>100504642   </td><td>Gm21996     </td><td>NA          </td></tr>
	<tr><th scope="row">27211</th><td>100041631   </td><td>Gm20879     </td><td>NA          </td></tr>
	<tr><th scope="row">27212</th><td>100504702   </td><td>NA          </td><td>NA          </td></tr>
	<tr><th scope="row">27213</th><td>100040357   </td><td>Gm20837     </td><td>NA          </td></tr>
	<tr><th scope="row">27214</th><td>100861808   </td><td>NA          </td><td>NA          </td></tr>
	<tr><th scope="row">27215</th><td>100504460   </td><td>NA          </td><td>NA          </td></tr>
	<tr><th scope="row">27216</th><td>100861837   </td><td>NA          </td><td>NA          </td></tr>
	<tr><th scope="row">27217</th><td>100861924   </td><td>NA          </td><td>NA          </td></tr>
	<tr><th scope="row">27218</th><td>170942      </td><td>Erdr1       </td><td>chrY        </td></tr>
	<tr><th scope="row">27219</th><td>100861691   </td><td>LOC100861691</td><td>NA          </td></tr>
	<tr><th scope="row">27220</th><td>100504472   </td><td>NA          </td><td>NA          </td></tr>
</tbody>
</table>
</dd>
</dl>



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
<font color ='#00bcd4'> In [11]: </font>

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
<font color ='#00bcd4'> In [12]: </font>

{% highlight R %}
table(rowSums(x$counts==0)==9)
{% endhighlight %}



    FALSE  TRUE
    22026  5153


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
<font color ='#00bcd4'> In [13]: </font>

{% highlight R %}
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

{% endhighlight %}


<ol class="list-inline">
	<li>14165</li>
	<li>9</li>
</ol>



Using this criterion, the number of genes is reduced to approximately half the
number that we started with (14,165 genes, panel B of the next figure). Note
that subsetting the entire DGEList-object removes both the counts as well as the
associated gene information. Code to produce the figure is given below.

<br>
<font color ='#00bcd4'> In [14]: </font>

{% highlight R %}
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")

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


![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_27_0.png)


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
<font color ='#00bcd4'> In [15]: </font>

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
<font color ='#00bcd4'> In [16]: </font>

{% highlight R %}
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
<font color ='#00bcd4'> In [17]: </font>

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




![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_34_1.png)


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
<font color ='#00bcd4'> In [18]: </font>

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


![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_38_0.png)


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
highlight cell population or sequencing lane (batch).

 An interactive MDS plot of this dataset can be found
<a href="{{site.url}}{{site.baseurl}}/MDSplot"> here </a>



<br>
<font color ='#00bcd4'> In [19]: </font>

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
<font color ='#00bcd4'> In [20]: </font>

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
<font color ='#00bcd4'> In [21]: </font>

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
<font color ='#00bcd4'> In [22]: </font>

{% highlight R %}
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
{% endhighlight %}


<dl>
	<dt>$genes</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">ENTREZID</th><th scope="col">SYMBOL</th><th scope="col">TXCHROM</th></tr></thead>
<tbody>
	<tr><th scope="row">1</th><td>497097 </td><td>Xkr4   </td><td>chr1   </td></tr>
	<tr><th scope="row">6</th><td>27395  </td><td>Mrpl15 </td><td>chr1   </td></tr>
	<tr><th scope="row">7</th><td>18777  </td><td>Lypla1 </td><td>chr1   </td></tr>
	<tr><th scope="row">9</th><td>21399  </td><td>Tcea1  </td><td>chr1   </td></tr>
	<tr><th scope="row">10</th><td>58175  </td><td>Rgs20  </td><td>chr1   </td></tr>
	<tr><th scope="row">11</th><td>108664 </td><td>Atp6v1h</td><td>chr1   </td></tr>
	<tr><th scope="row">14</th><td>12421  </td><td>Rb1cc1 </td><td>chr1   </td></tr>
	<tr><th scope="row">17</th><td>319263 </td><td>Pcmtd1 </td><td>chr1   </td></tr>
	<tr><th scope="row">19</th><td>59014  </td><td>Rrs1   </td><td>chr1   </td></tr>
	<tr><th scope="row">20</th><td>76187  </td><td>Adhfe1 </td><td>chr1   </td></tr>
	<tr><th scope="row">23</th><td>17864  </td><td>Mybl1  </td><td>chr1   </td></tr>
	<tr><th scope="row">24</th><td>70675  </td><td>Vcpip1 </td><td>chr1   </td></tr>
	<tr><th scope="row">26</th><td>170755 </td><td>Sgk3   </td><td>chr1   </td></tr>
	<tr><th scope="row">27</th><td>620986 </td><td>Gm6195 </td><td>NA     </td></tr>
	<tr><th scope="row">29</th><td>73824  </td><td>Snhg6  </td><td>chr1   </td></tr>
	<tr><th scope="row">34</th><td>26754  </td><td>Cops5  </td><td>chr1   </td></tr>
	<tr><th scope="row">35</th><td>211660 </td><td>Cspp1  </td><td>chr1   </td></tr>
	<tr><th scope="row">36</th><td>211673 </td><td>Arfgef1</td><td>chr1   </td></tr>
	<tr><th scope="row">37</th><td>329093 </td><td>Cpa6   </td><td>chr1   </td></tr>
	<tr><th scope="row">38</th><td>109294 </td><td>Prex2  </td><td>chr1   </td></tr>
	<tr><th scope="row">41</th><td>240725 </td><td>Sulf1  </td><td>chr1   </td></tr>
	<tr><th scope="row">44</th><td>17978  </td><td>Ncoa2  </td><td>chr1   </td></tr>
	<tr><th scope="row">45</th><td>72265  </td><td>Tram1  </td><td>chr1   </td></tr>
	<tr><th scope="row">46</th><td>212442 </td><td>Lactb2 </td><td>chr1   </td></tr>
	<tr><th scope="row">49</th><td>14048  </td><td>Eya1   </td><td>chr1   </td></tr>
	<tr><th scope="row">50</th><td>17681  </td><td>Msc    </td><td>chr1   </td></tr>
	<tr><th scope="row">53</th><td>21749  </td><td>Terf1  </td><td>chr1   </td></tr>
	<tr><th scope="row">54</th><td>226866 </td><td>Sbspon </td><td>chr1   </td></tr>
	<tr><th scope="row">57</th><td>19989  </td><td>Rpl7   </td><td>chr1   </td></tr>
	<tr><th scope="row">58</th><td>98711  </td><td>Rdh10  </td><td>chr1   </td></tr>
	<tr><th scope="row">⋮</th><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope="row">26853</th><td>245688   </td><td>Rbbp7    </td><td>chrX     </td></tr>
	<tr><th scope="row">26854</th><td>353170   </td><td>Txlng    </td><td>chrX     </td></tr>
	<tr><th scope="row">26855</th><td>67043    </td><td>Syap1    </td><td>chrX     </td></tr>
	<tr><th scope="row">26858</th><td>55936    </td><td>Ctps2    </td><td>chrX     </td></tr>
	<tr><th scope="row">26862</th><td>108012   </td><td>Ap1s2    </td><td>chrX     </td></tr>
	<tr><th scope="row">26863</th><td>22184    </td><td>Zrsr2    </td><td>chrX     </td></tr>
	<tr><th scope="row">26864</th><td>56078    </td><td>Car5b    </td><td>chrX     </td></tr>
	<tr><th scope="row">26865</th><td>20438    </td><td>Siah1b   </td><td>chrX     </td></tr>
	<tr><th scope="row">26868</th><td>12169    </td><td>Bmx      </td><td>chrX     </td></tr>
	<tr><th scope="row">26869</th><td>69656    </td><td>Pir      </td><td>chrX     </td></tr>
	<tr><th scope="row">26870</th><td>14205    </td><td>Vegfd    </td><td>chrX     </td></tr>
	<tr><th scope="row">26871</th><td>18700    </td><td>Piga     </td><td>chrX     </td></tr>
	<tr><th scope="row">26874</th><td>76763    </td><td>Mospd2   </td><td>chrX     </td></tr>
	<tr><th scope="row">26875</th><td>237211   </td><td>Fancb    </td><td>chrX     </td></tr>
	<tr><th scope="row">26877</th><td>237221   </td><td>Gemin8   </td><td>chrX     </td></tr>
	<tr><th scope="row">26878</th><td>14758    </td><td>Gpm6b    </td><td>chrX     </td></tr>
	<tr><th scope="row">26879</th><td>237222   </td><td>Ofd1     </td><td>chrX     </td></tr>
	<tr><th scope="row">26880</th><td>66226    </td><td>Trappc2  </td><td>chrX     </td></tr>
	<tr><th scope="row">26881</th><td>56382    </td><td>Rab9     </td><td>chrX     </td></tr>
	<tr><th scope="row">26882</th><td>245695   </td><td>Tceanc   </td><td>chrX     </td></tr>
	<tr><th scope="row">26883</th><td>54156    </td><td>Egfl6    </td><td>chrX     </td></tr>
	<tr><th scope="row">26891</th><td>110639   </td><td>Prps2    </td><td>chrX     </td></tr>
	<tr><th scope="row">26893</th><td>333605   </td><td>Frmpd4   </td><td>chrX     </td></tr>
	<tr><th scope="row">26894</th><td>17692    </td><td>Msl3     </td><td>chrX     </td></tr>
	<tr><th scope="row">26895</th><td>11856    </td><td>Arhgap6  </td><td>chrX     </td></tr>
	<tr><th scope="row">26897</th><td>15159    </td><td>Hccs     </td><td>chrX     </td></tr>
	<tr><th scope="row">26899</th><td>17318    </td><td>Mid1     </td><td>chrX     </td></tr>
	<tr><th scope="row">26903</th><td>100862004</td><td>NA       </td><td>NA       </td></tr>
	<tr><th scope="row">27216</th><td>100861837</td><td>NA       </td><td>NA       </td></tr>
	<tr><th scope="row">27218</th><td>170942   </td><td>Erdr1    </td><td>chrY     </td></tr>
</tbody>
</table>
</dd>
	<dt>$targets</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">files</th><th scope="col">group</th><th scope="col">lib.size</th><th scope="col">norm.factors</th><th scope="col">lane</th></tr></thead>
<tbody>
	<tr><th scope="row">10_6_5_11</th><td>GSM1545535_10_6_5_11.txt</td><td>LP                      </td><td>29409426                </td><td>0.8957309               </td><td>L004                    </td></tr>
	<tr><th scope="row">9_6_5_11</th><td>GSM1545536_9_6_5_11.txt </td><td>ML                      </td><td>36528591                </td><td>1.0349196               </td><td>L004                    </td></tr>
	<tr><th scope="row">purep53</th><td>GSM1545538_purep53.txt  </td><td>Basal                   </td><td>59598629                </td><td>1.0439552               </td><td>L004                    </td></tr>
	<tr><th scope="row">JMS8-2</th><td>GSM1545539_JMS8-2.txt   </td><td>Basal                   </td><td>53382070                </td><td>1.0405040               </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-3</th><td>GSM1545540_JMS8-3.txt   </td><td>ML                      </td><td>78175314                </td><td>1.0323599               </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-4</th><td>GSM1545541_JMS8-4.txt   </td><td>LP                      </td><td>55762781                </td><td>0.9223424               </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS8-5</th><td>GSM1545542_JMS8-5.txt   </td><td>Basal                   </td><td>54115150                </td><td>0.9836603               </td><td>L006                    </td></tr>
	<tr><th scope="row">JMS9-P7c</th><td>GSM1545544_JMS9-P7c.txt </td><td>ML                      </td><td>23043111                </td><td>1.0827381               </td><td>L008                    </td></tr>
	<tr><th scope="row">JMS9-P8c</th><td>GSM1545545_JMS9-P8c.txt </td><td>LP                      </td><td>19525423                </td><td>0.9792607               </td><td>L008                    </td></tr>
</tbody>
</table>
</dd>
	<dt>$E</dt>
		<dd><table class="table-responsive table-striped">
<thead><tr><th></th><th scope="col">10_6_5_11</th><th scope="col">9_6_5_11</th><th scope="col">purep53</th><th scope="col">JMS8-2</th><th scope="col">JMS8-3</th><th scope="col">JMS8-4</th><th scope="col">JMS8-5</th><th scope="col">JMS9-P7c</th><th scope="col">JMS9-P8c</th></tr></thead>
<tbody>
	<tr><th scope="row">497097</th><td>-4.2932443 </td><td>-3.869026  </td><td> 2.5227529 </td><td> 3.3020063 </td><td>-4.4812863 </td><td>-3.99387571</td><td> 3.3067821 </td><td>-3.2043356 </td><td>-5.2872819 </td></tr>
	<tr><th scope="row">27395</th><td> 3.8750100 </td><td> 4.400568  </td><td> 4.5211725 </td><td> 4.5706244 </td><td> 4.3228447 </td><td> 3.78654688</td><td> 3.9188779 </td><td> 4.3456416 </td><td> 4.1326782 </td></tr>
	<tr><th scope="row">18777</th><td> 4.7076947 </td><td> 5.559334  </td><td> 5.4005688 </td><td> 5.1712347 </td><td> 5.6277981 </td><td> 5.08179444</td><td> 5.0800614 </td><td> 5.7574035 </td><td> 5.1504701 </td></tr>
	<tr><th scope="row">21399</th><td> 4.7844616 </td><td> 4.741999  </td><td> 5.3745475 </td><td> 5.1309249 </td><td> 4.8480295 </td><td> 4.94402351</td><td> 5.1582920 </td><td> 5.0369326 </td><td> 4.9876785 </td></tr>
	<tr><th scope="row">58175</th><td> 3.9435672 </td><td> 3.294875  </td><td>-1.7679243 </td><td>-1.8803024 </td><td> 2.9932888 </td><td> 3.35737906</td><td>-2.1141045 </td><td> 3.1426213 </td><td> 3.5232897 </td></tr>
	<tr><th scope="row">108664</th><td> 5.8670474 </td><td> 6.196255  </td><td> 5.1372478 </td><td> 5.2781767 </td><td> 6.1364435 </td><td> 6.02472190</td><td> 4.9998460 </td><td> 6.2181498 </td><td> 6.0021606 </td></tr>
	<tr><th scope="row">12421</th><td> 6.8748010 </td><td> 6.207522  </td><td> 5.8338996 </td><td> 5.6692500 </td><td> 6.0649192 </td><td> 6.81543307</td><td> 6.0747316 </td><td> 6.3398561 </td><td> 6.7642669 </td></tr>
	<tr><th scope="row">319263</th><td> 6.1065677 </td><td> 5.798795  </td><td> 6.0396147 </td><td> 5.6687188 </td><td> 5.6307808 </td><td> 6.04758789</td><td> 6.0532133 </td><td> 6.1261332 </td><td> 6.1098579 </td></tr>
	<tr><th scope="row">59014</th><td> 5.0234144 </td><td> 4.753758  </td><td> 6.5273527 </td><td> 6.3028624 </td><td> 5.4689571 </td><td> 5.71717664</td><td> 6.6237175 </td><td> 4.8169221 </td><td> 5.2274321 </td></tr>
	<tr><th scope="row">76187</th><td> 0.8899776 </td><td> 1.066434  </td><td>-0.3580485 </td><td>-0.7159156 </td><td> 0.9256779 </td><td> 1.18176294</td><td> 0.1368571 </td><td> 0.9494698 </td><td> 0.4940778 </td></tr>
	<tr><th scope="row">17864</th><td> 2.7254196 </td><td> 2.188424  </td><td> 4.0793569 </td><td> 4.0140973 </td><td> 3.0162804 </td><td> 3.41430237</td><td> 4.0582230 </td><td> 2.1247880 </td><td> 2.7516370 </td></tr>
	<tr><th scope="row">70675</th><td> 6.3448899 </td><td> 5.783819  </td><td> 5.8761380 </td><td> 5.8222883 </td><td> 6.1134382 </td><td> 6.30353191</td><td> 6.2690814 </td><td> 6.1044585 </td><td> 6.2849446 </td></tr>
	<tr><th scope="row">170755</th><td> 1.5730043 </td><td> 3.146668  </td><td> 3.6234114 </td><td> 3.4050998 </td><td> 3.2062144 </td><td> 2.16166537</td><td> 3.6225004 </td><td> 3.5371314 </td><td> 2.5392665 </td></tr>
	<tr><th scope="row">620986</th><td> 1.2816646 </td><td> 2.412672  </td><td> 2.6858755 </td><td> 1.9517146 </td><td> 1.9758014 </td><td> 0.98667193</td><td> 1.9320373 </td><td> 2.2875175 </td><td> 1.3709295 </td></tr>
	<tr><th scope="row">73824</th><td> 3.2992128 </td><td> 3.738304  </td><td> 3.7901684 </td><td> 3.4050998 </td><td> 3.1820176 </td><td> 3.21418442</td><td> 3.6138160 </td><td> 3.8639053 </td><td> 3.7975264 </td></tr>
	<tr><th scope="row">26754</th><td> 5.9391767 </td><td> 6.096470  </td><td> 6.0948459 </td><td> 5.9970611 </td><td> 6.1484606 </td><td> 5.90634420</td><td> 5.8026110 </td><td> 6.2315430 </td><td> 6.0170694 </td></tr>
	<tr><th scope="row">211660</th><td> 4.6638578 </td><td> 4.499917  </td><td> 4.5264339 </td><td> 4.2397120 </td><td> 4.4766449 </td><td> 5.26047807</td><td> 4.9695344 </td><td> 5.0426425 </td><td> 4.9992758 </td></tr>
	<tr><th scope="row">211673</th><td> 7.8309848 </td><td> 7.609036  </td><td> 7.5756739 </td><td> 7.2225331 </td><td> 7.7648180 </td><td> 8.08375070</td><td> 7.7237757 </td><td> 7.7488424 </td><td> 7.7695251 </td></tr>
	<tr><th scope="row">329093</th><td>-3.5562787 </td><td>-6.190954  </td><td> 3.7511503 </td><td> 3.3741561 </td><td>-7.2886412 </td><td>-3.99387571</td><td> 3.5099964 </td><td>-5.5262637 </td><td>-3.7023194 </td></tr>
	<tr><th scope="row">109294</th><td> 0.4968327 </td><td> 1.371288  </td><td> 2.4843357 </td><td> 2.1716097 </td><td> 2.2368796 </td><td> 2.39106218</td><td> 3.6839460 </td><td> 2.2084460 </td><td> 1.8107501 </td></tr>
	<tr><th scope="row">240725</th><td> 1.8565029 </td><td> 1.502533  </td><td> 5.2155580 </td><td> 5.2336194 </td><td> 0.9640242 </td><td> 2.05985627</td><td> 5.7915834 </td><td> 1.0736492 </td><td> 1.4269636 </td></tr>
	<tr><th scope="row">17978</th><td> 7.2399232 </td><td> 7.274995  </td><td> 6.8600786 </td><td> 6.8961874 </td><td> 7.3372388 </td><td> 7.15868077</td><td> 7.0915422 </td><td> 7.4346433 </td><td> 7.2501792 </td></tr>
	<tr><th scope="row">72265</th><td> 7.5318470 </td><td> 7.200155  </td><td> 6.7060713 </td><td> 6.7408653 </td><td> 7.4346464 </td><td> 7.93965688</td><td> 6.8273540 </td><td> 7.0155590 </td><td> 7.5701160 </td></tr>
	<tr><th scope="row">212442</th><td> 3.5076556 </td><td> 4.607518  </td><td> 2.6707488 </td><td> 2.7313584 </td><td> 4.3619618 </td><td> 3.92754022</td><td> 2.3583833 </td><td> 4.9645872 </td><td> 4.0278676 </td></tr>
	<tr><th scope="row">14048</th><td>-4.2932443 </td><td>-4.605992  </td><td> 4.2569779 </td><td> 4.2296633 </td><td>-3.2011784 </td><td>-1.18652079</td><td> 4.1315433 </td><td>-1.2783361 </td><td>-0.8949645 </td></tr>
	<tr><th scope="row">17681</th><td>-0.1502863 </td><td> 2.613177  </td><td>-2.9903167 </td><td>-4.4163553 </td><td> 1.4290352 </td><td>-0.46138063</td><td>-2.6704978 </td><td> 2.2875175 </td><td> 0.6899980 </td></tr>
	<tr><th scope="row">21749</th><td> 3.1797850 </td><td> 3.361715  </td><td> 3.7866646 </td><td> 3.2854710 </td><td> 3.0185596 </td><td> 2.99805099</td><td> 3.3201902 </td><td> 3.9835113 </td><td> 3.4775896 </td></tr>
	<tr><th scope="row">226866</th><td>-1.7907439 </td><td>-0.905552  </td><td> 6.2498386 </td><td> 6.1456963 </td><td>-1.1388941 </td><td>-0.05976365</td><td> 6.1927775 </td><td>-0.5720673 </td><td>-0.8949645 </td></tr>
	<tr><th scope="row">19989</th><td> 9.2243021 </td><td> 8.817518  </td><td> 9.2143156 </td><td> 8.8301775 </td><td> 8.6309437 </td><td> 8.76497209</td><td> 8.8795982 </td><td> 8.3692170 </td><td> 8.7607791 </td></tr>
	<tr><th scope="row">98711</th><td> 5.6864191 </td><td> 5.457853  </td><td> 6.0259345 </td><td> 6.3225814 </td><td> 5.1373614 </td><td> 5.65948122</td><td> 5.6140882 </td><td> 4.9625807 </td><td> 5.0492246 </td></tr>
	<tr><th scope="row">⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope="row">245688</th><td> 6.6880852  </td><td> 6.6393629  </td><td> 8.0241223  </td><td> 7.9215453  </td><td> 7.153589557</td><td> 6.9468574  </td><td>7.81254734  </td><td> 6.58389336 </td><td> 6.4993961  </td></tr>
	<tr><th scope="row">353170</th><td> 4.0029072  </td><td> 3.3306462  </td><td> 3.1862721  </td><td> 2.8143857  </td><td> 3.470414714</td><td> 4.1924154  </td><td>3.59849131  </td><td> 3.97954789 </td><td> 4.1157301  </td></tr>
	<tr><th scope="row">67043</th><td> 6.2162009  </td><td> 6.4852232  </td><td> 6.0340843  </td><td> 5.6585889  </td><td> 6.423993243</td><td> 6.5185826  </td><td>6.26527390  </td><td> 6.74198591 </td><td> 6.4879168  </td></tr>
	<tr><th scope="row">55936</th><td> 6.1678941  </td><td> 6.0767103  </td><td> 5.1550211  </td><td> 5.3119054  </td><td> 6.132240088</td><td> 6.1054721  </td><td>5.13950668  </td><td> 6.10445850 </td><td> 6.1363593  </td></tr>
	<tr><th scope="row">108012</th><td> 1.9483417  </td><td> 0.3789014  </td><td> 3.7257597  </td><td> 3.1397675  </td><td> 1.019697805</td><td> 2.5184415  </td><td>4.35121705  </td><td> 0.40447368 </td><td> 1.2825737  </td></tr>
	<tr><th scope="row">22184</th><td> 5.7323567  </td><td> 5.1076812  </td><td> 5.6113303  </td><td> 5.0702787  </td><td> 5.133686399</td><td> 5.5377842  </td><td>5.55520582  </td><td> 5.46313084 </td><td> 5.5240928  </td></tr>
	<tr><th scope="row">56078</th><td> 4.4175622  </td><td> 2.5267222  </td><td> 1.7787498  </td><td> 2.0722883  </td><td> 3.436725033</td><td> 4.8212788  </td><td>2.00691093  </td><td> 2.30028483 </td><td> 4.2025660  </td></tr>
	<tr><th scope="row">20438</th><td> 2.6096333  </td><td> 2.6639141  </td><td> 3.4031453  </td><td> 3.2231663  </td><td> 2.658265050</td><td> 2.6060371  </td><td>2.94421203  </td><td> 2.57702415 </td><td> 2.2440995  </td></tr>
	<tr><th scope="row">12169</th><td>-4.2932443  </td><td>-1.4360668  </td><td>-0.2389958  </td><td> 0.1072067  </td><td>-1.138894105</td><td>-1.1865208  </td><td>0.08752939  </td><td> 0.02832519 </td><td>-1.1998191  </td></tr>
	<tr><th scope="row">69656</th><td> 2.1607122  </td><td> 4.7654219  </td><td> 3.3447759  </td><td> 2.8634874  </td><td> 3.951553829</td><td> 3.3395991  </td><td>2.78700378  </td><td> 5.00998356 </td><td> 3.8901376  </td></tr>
	<tr><th scope="row">14205</th><td>-1.1233193  </td><td> 0.5232913  </td><td>-0.7073827  </td><td>-0.2953399  </td><td> 0.005979524</td><td> 0.9334790  </td><td>0.27546234  </td><td> 1.39259958 </td><td> 0.4940778  </td></tr>
	<tr><th scope="row">18700</th><td> 3.7235640  </td><td> 4.0770028  </td><td> 4.4420860  </td><td> 4.0746961  </td><td> 4.025942019</td><td> 4.0882733  </td><td>4.61327147  </td><td> 4.72521875 </td><td> 4.4218019  </td></tr>
	<tr><th scope="row">76763</th><td> 5.2822950  </td><td> 5.0961809  </td><td> 3.8915111  </td><td> 3.5620692  </td><td> 4.884098792</td><td> 5.3777454  </td><td>4.06620255  </td><td> 5.51334086 </td><td> 5.5589920  </td></tr>
	<tr><th scope="row">237211</th><td> 0.8899776  </td><td> 0.6292247  </td><td> 2.0006382  </td><td> 1.3115652  </td><td> 0.273601199</td><td> 0.3786785  </td><td>1.09378838  </td><td> 0.91667984 </td><td> 1.1556616  </td></tr>
	<tr><th scope="row">237221</th><td> 2.0933368  </td><td> 1.6228269  </td><td> 2.1876011  </td><td> 2.4490687  </td><td> 1.462902834</td><td> 0.9866719  </td><td>1.78507116  </td><td> 1.90836457 </td><td> 2.0436349  </td></tr>
	<tr><th scope="row">14758</th><td> 5.9343717  </td><td> 4.2509525  </td><td> 5.2978573  </td><td> 5.2493359  </td><td> 5.117829544</td><td> 6.0452392  </td><td>5.42163747  </td><td> 3.87674836 </td><td> 5.7399329  </td></tr>
	<tr><th scope="row">237222</th><td> 4.4605296  </td><td> 4.1093983  </td><td> 4.2595079  </td><td> 4.0307279  </td><td> 4.575157879</td><td> 4.7719432  </td><td>4.55435509  </td><td> 4.52494528 </td><td> 4.6508274  </td></tr>
	<tr><th scope="row">66226</th><td> 0.7511499  </td><td> 1.1219287  </td><td> 1.0743363  </td><td> 0.6454209  </td><td> 1.133423541</td><td> 1.1931228  </td><td>1.15492868  </td><td> 1.41625085 </td><td> 0.9415367  </td></tr>
	<tr><th scope="row">56382</th><td> 4.8386127  </td><td> 4.2776698  </td><td> 4.3453714  </td><td> 4.2666393  </td><td> 4.244201686</td><td> 4.5133526  </td><td>4.17941316  </td><td> 4.72284979 </td><td> 4.7557453  </td></tr>
	<tr><th scope="row">245695</th><td> 1.7728449  </td><td> 2.5936806  </td><td> 0.6184926  </td><td> 0.8391454  </td><td> 1.615240621</td><td> 1.2594653  </td><td>0.97674896  </td><td> 2.37460315 </td><td> 1.6075358  </td></tr>
	<tr><th scope="row">54156</th><td>-0.7489237  </td><td>-1.5470981  </td><td> 5.5813096  </td><td> 5.3294871  </td><td>-0.574395707</td><td> 1.8023957  </td><td>5.43031910  </td><td>-1.27833615 </td><td> 0.6899980  </td></tr>
	<tr><th scope="row">110639</th><td> 5.6300832  </td><td> 5.1846283  </td><td> 6.1771015  </td><td> 6.1959604  </td><td> 5.619563728</td><td> 5.6486918  </td><td>6.02606067  </td><td> 5.43590937 </td><td> 5.5853929  </td></tr>
	<tr><th scope="row">333605</th><td>-5.8782068  </td><td>-3.3835993  </td><td> 1.7716777  </td><td> 1.7969920  </td><td>-5.703678724</td><td>-4.4793025  </td><td>1.38669759  </td><td>-2.71890874 </td><td>-3.7023194  </td></tr>
	<tr><th scope="row">17692</th><td> 5.3355050  </td><td> 4.9883329  </td><td> 5.3533872  </td><td> 5.6239344  </td><td> 5.198947622</td><td> 5.1885186  </td><td>5.24309562  </td><td> 5.08199116 </td><td> 5.2215032  </td></tr>
	<tr><th scope="row">11856</th><td> 5.6739820  </td><td> 3.9214853  </td><td> 3.5153626  </td><td> 3.1850441  </td><td> 4.487380490</td><td> 5.2290917  </td><td>3.32551867  </td><td> 3.74753194 </td><td> 5.5869310  </td></tr>
	<tr><th scope="row">15159</th><td> 4.0421461  </td><td> 3.5555601  </td><td> 4.7416811  </td><td> 4.2382807  </td><td> 4.130792326</td><td> 4.8028592  </td><td>4.66357825  </td><td> 3.94744209 </td><td> 4.4044616  </td></tr>
	<tr><th scope="row">17318</th><td> 4.1232014  </td><td> 5.4968588  </td><td> 4.1020881  </td><td> 3.8846836  </td><td> 5.026225206</td><td> 3.9825861  </td><td>4.38733468  </td><td> 5.54586894 </td><td> 4.0050397  </td></tr>
	<tr><th scope="row">100862004</th><td> 5.4454110  </td><td> 5.1202264  </td><td> 5.5166851  </td><td> 4.8722801  </td><td> 4.971396196</td><td> 5.1504184  </td><td>5.51903608  </td><td> 5.09120381 </td><td> 4.9900055  </td></tr>
	<tr><th scope="row">100861837</th><td> 3.4459738  </td><td> 3.6932163  </td><td> 3.2889070  </td><td> 2.6518856  </td><td> 3.375806060</td><td> 2.9553257  </td><td>3.20926560  </td><td> 3.12119477 </td><td> 2.1135975  </td></tr>
	<tr><th scope="row">170942</th><td> 4.8091689  </td><td> 4.3968233  </td><td> 4.4881159  </td><td> 3.8718187  </td><td> 4.017990137</td><td> 4.2882198  </td><td>4.55775500  </td><td> 4.34871769 </td><td> 4.2806741  </td></tr>
</tbody>
</table>
</dd>
	<dt>$weights</dt>
		<dd><table class="table-responsive table-striped">
<tbody>
	<tr><td> 1.183974</td><td> 1.183974</td><td>20.526779</td><td>20.977471</td><td> 1.773562</td><td> 1.217142</td><td>21.125740</td><td> 1.183974</td><td> 1.183974</td></tr>
	<tr><td>20.879554</td><td>26.561871</td><td>31.596323</td><td>29.661022</td><td>32.558344</td><td>26.745293</td><td>29.792090</td><td>21.900102</td><td>17.150677</td></tr>
	<tr><td>28.003202</td><td>33.695540</td><td>34.845507</td><td>34.456731</td><td>35.148529</td><td>33.550527</td><td>34.517259</td><td>31.440457</td><td>25.228325</td></tr>
	<tr><td>27.670233</td><td>29.595778</td><td>34.901302</td><td>34.432980</td><td>34.841349</td><td>33.159425</td><td>34.493456</td><td>26.136796</td><td>24.502247</td></tr>
	<tr><td>19.737381</td><td>18.658333</td><td> 3.184207</td><td> 2.629860</td><td>24.191635</td><td>24.014937</td><td> 2.648747</td><td>13.149278</td><td>14.351930</td></tr>
	<tr><td>33.902057</td><td>35.238072</td><td>34.618260</td><td>34.210465</td><td>34.263943</td><td>35.239031</td><td>34.282857</td><td>33.648502</td><td>31.287895</td></tr>
	<tr><td>35.273872</td><td>35.296976</td><td>35.220307</td><td>35.325989</td><td>34.309998</td><td>34.014915</td><td>35.317151</td><td>33.748949</td><td>34.941572</td></tr>
	<tr><td>34.453035</td><td>34.689522</td><td>35.123118</td><td>35.312597</td><td>35.027957</td><td>35.252747</td><td>35.303763</td><td>32.606997</td><td>32.591533</td></tr>
	<tr><td>30.004066</td><td>30.023974</td><td>34.883075</td><td>34.498348</td><td>35.328587</td><td>35.294917</td><td>34.461700</td><td>25.126350</td><td>25.668016</td></tr>
	<tr><td> 6.020713</td><td> 7.298922</td><td> 5.386546</td><td> 5.193771</td><td>12.050396</td><td> 9.336468</td><td> 5.236998</td><td> 4.899833</td><td> 4.225656</td></tr>
	<tr><td>13.822179</td><td>12.752437</td><td>27.126863</td><td>29.090460</td><td>22.708735</td><td>22.963293</td><td>29.218382</td><td> 8.989551</td><td>10.153537</td></tr>
	<tr><td>34.832186</td><td>34.836659</td><td>35.222840</td><td>35.226472</td><td>34.513910</td><td>34.807476</td><td>35.211027</td><td>32.816510</td><td>33.148596</td></tr>
	<tr><td> 9.059951</td><td>16.672560</td><td>24.761013</td><td>24.694918</td><td>25.858365</td><td>14.038225</td><td>24.851859</td><td>16.211528</td><td> 9.062693</td></tr>
	<tr><td> 7.501635</td><td>13.112928</td><td>17.939276</td><td>13.822286</td><td>16.421315</td><td> 8.976161</td><td>13.928057</td><td> 9.628375</td><td> 5.559286</td></tr>
	<tr><td>17.125103</td><td>20.418795</td><td>27.486704</td><td>23.923279</td><td>26.516004</td><td>21.524498</td><td>24.074110</td><td>17.966622</td><td>15.428681</td></tr>
	<tr><td>34.001221</td><td>35.249141</td><td>35.123544</td><td>35.282639</td><td>34.457084</td><td>35.279785</td><td>35.273815</td><td>33.676076</td><td>31.411303</td></tr>
	<tr><td>27.101197</td><td>27.186150</td><td>31.591640</td><td>31.852162</td><td>34.461645</td><td>33.778849</td><td>31.960279</td><td>25.147432</td><td>25.627445</td></tr>
	<tr><td>33.558537</td><td>33.218843</td><td>31.862269</td><td>32.099482</td><td>29.546789</td><td>30.562329</td><td>32.038941</td><td>34.776429</td><td>34.850318</td></tr>
	<tr><td> 1.183974</td><td> 1.183974</td><td>28.390218</td><td>22.749409</td><td> 1.183974</td><td> 1.267314</td><td>22.900138</td><td> 1.183974</td><td> 1.183974</td></tr>
	<tr><td> 6.068260</td><td> 8.212458</td><td>15.237872</td><td>20.964463</td><td>19.122873</td><td>13.746616</td><td>21.113229</td><td> 9.060281</td><td> 6.913503</td></tr>
	<tr><td> 9.232370</td><td> 8.133480</td><td>35.192808</td><td>34.965256</td><td>12.976879</td><td>13.562400</td><td>35.005537</td><td> 5.250756</td><td> 6.134439</td></tr>
	<tr><td>34.978931</td><td>34.224620</td><td>33.392427</td><td>33.687323</td><td>30.976639</td><td>32.851825</td><td>33.636238</td><td>35.097007</td><td>35.327319</td></tr>
	<tr><td>34.134629</td><td>34.464042</td><td>34.178124</td><td>34.027061</td><td>30.766592</td><td>30.844409</td><td>33.981483</td><td>35.323848</td><td>35.200487</td></tr>
	<tr><td>18.876545</td><td>27.483863</td><td>18.290361</td><td>17.228688</td><td>33.873418</td><td>25.882132</td><td>17.359868</td><td>25.436728</td><td>17.534849</td></tr>
	<tr><td> 1.191804</td><td> 1.183974</td><td>23.783125</td><td>31.856364</td><td> 2.182204</td><td> 2.531202</td><td>31.964007</td><td> 2.042838</td><td> 2.586535</td></tr>
	<tr><td> 4.367502</td><td>12.154050</td><td> 2.040621</td><td> 1.542420</td><td>14.356298</td><td> 4.842211</td><td> 1.552412</td><td>10.461685</td><td> 3.885675</td></tr>
	<tr><td>15.644306</td><td>19.282242</td><td>26.675346</td><td>22.724600</td><td>24.949989</td><td>19.448771</td><td>22.875125</td><td>17.556056</td><td>14.616158</td></tr>
	<tr><td> 2.607327</td><td> 2.971384</td><td>35.197913</td><td>34.942430</td><td> 5.406177</td><td> 4.342430</td><td>34.917508</td><td> 2.743770</td><td> 2.476632</td></tr>
	<tr><td>29.533758</td><td>29.527115</td><td>26.069536</td><td>27.697187</td><td>26.952972</td><td>27.539516</td><td>27.630875</td><td>33.064485</td><td>32.862533</td></tr>
	<tr><td>32.821051</td><td>32.976639</td><td>35.058318</td><td>35.280685</td><td>35.328263</td><td>35.260148</td><td>35.271861</td><td>25.118022</td><td>25.424523</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>35.313518</td><td>35.020062</td><td>30.738346</td><td>30.705088</td><td>32.106349</td><td>33.690363</td><td>30.640712</td><td>34.863909</td><td>33.817288</td></tr>
	<tr><td>20.779410</td><td>19.336050</td><td>22.270921</td><td>21.867296</td><td>28.420824</td><td>28.583229</td><td>22.018730</td><td>17.182757</td><td>19.040610</td></tr>
	<tr><td>35.043275</td><td>35.294935</td><td>35.167400</td><td>35.255865</td><td>33.451692</td><td>34.818201</td><td>35.242531</td><td>35.077315</td><td>33.985950</td></tr>
	<tr><td>34.575479</td><td>35.175701</td><td>34.810060</td><td>34.413720</td><td>34.403484</td><td>35.079156</td><td>34.474154</td><td>33.117615</td><td>32.062429</td></tr>
	<tr><td> 9.442684</td><td> 6.070304</td><td>25.088249</td><td>26.716679</td><td>11.482491</td><td>15.861477</td><td>26.862484</td><td> 3.852298</td><td> 6.003273</td></tr>
	<tr><td>32.588210</td><td>32.318592</td><td>35.309374</td><td>34.834493</td><td>35.299755</td><td>35.168698</td><td>34.874586</td><td>28.436238</td><td>29.271128</td></tr>
	<tr><td>24.330267</td><td>14.501332</td><td>12.604272</td><td>14.377398</td><td>25.756255</td><td>33.509253</td><td>14.486847</td><td> 9.848796</td><td>17.734866</td></tr>
	<tr><td>12.540541</td><td>14.984426</td><td>23.415665</td><td>21.200781</td><td>21.544249</td><td>17.040313</td><td>21.346462</td><td>10.259992</td><td> 8.680450</td></tr>
	<tr><td> 1.288990</td><td> 2.269900</td><td> 4.516606</td><td> 6.908047</td><td> 5.398646</td><td> 2.609343</td><td> 6.968533</td><td> 3.730206</td><td> 2.060616</td></tr>
	<tr><td>13.839189</td><td>26.084909</td><td>21.772162</td><td>19.862781</td><td>32.555236</td><td>18.995981</td><td>20.003657</td><td>27.453991</td><td>15.283409</td></tr>
	<tr><td> 3.316814</td><td> 4.633936</td><td> 4.607345</td><td> 5.885335</td><td>10.274683</td><td> 6.547269</td><td> 5.935741</td><td> 5.876179</td><td> 4.250339</td></tr>
	<tr><td>20.052442</td><td>23.935182</td><td>30.997310</td><td>30.358081</td><td>31.905705</td><td>27.524314</td><td>30.479328</td><td>23.609976</td><td>20.279926</td></tr>
	<tr><td>30.722861</td><td>31.144850</td><td>28.032769</td><td>26.389211</td><td>35.188842</td><td>34.722473</td><td>26.532783</td><td>29.134402</td><td>29.152740</td></tr>
	<tr><td> 6.063701</td><td> 6.364891</td><td>13.802363</td><td>10.276149</td><td> 8.004177</td><td> 7.052066</td><td>10.360417</td><td> 5.218934</td><td> 5.137554</td></tr>
	<tr><td> 8.760525</td><td> 9.849421</td><td>17.123929</td><td>13.753105</td><td>13.034399</td><td>10.923870</td><td>13.858214</td><td> 8.146730</td><td> 7.468315</td></tr>
	<tr><td>33.890024</td><td>26.919064</td><td>34.615106</td><td>34.977651</td><td>34.676024</td><td>34.981691</td><td>35.017950</td><td>19.033151</td><td>28.457889</td></tr>
	<tr><td>24.756618</td><td>25.439840</td><td>29.613203</td><td>30.166352</td><td>33.778657</td><td>32.552068</td><td>30.286643</td><td>22.084990</td><td>21.973762</td></tr>
	<tr><td> 5.939548</td><td> 7.649549</td><td> 9.419880</td><td> 9.046286</td><td>12.530571</td><td> 9.184240</td><td> 9.124145</td><td> 6.129156</td><td> 4.949981</td></tr>
	<tr><td>26.480629</td><td>26.664618</td><td>31.157965</td><td>29.188882</td><td>32.654208</td><td>31.562417</td><td>29.317345</td><td>23.003548</td><td>23.383993</td></tr>
	<tr><td> 8.327938</td><td>12.473099</td><td> 9.912593</td><td> 7.997338</td><td>16.779984</td><td>10.748803</td><td> 8.067718</td><td> 9.545503</td><td> 6.454228</td></tr>
	<tr><td> 4.282901</td><td> 2.518129</td><td>33.410863</td><td>35.320617</td><td> 5.661246</td><td> 9.702991</td><td>35.323130</td><td> 2.340881</td><td> 4.067437</td></tr>
	<tr><td>32.191084</td><td>32.750595</td><td>35.072219</td><td>35.125383</td><td>35.241158</td><td>35.316777</td><td>35.105387</td><td>29.078450</td><td>28.878699</td></tr>
	<tr><td> 1.183974</td><td> 1.183974</td><td>13.959083</td><td>11.441190</td><td> 1.301688</td><td> 1.183974</td><td>11.532998</td><td> 1.327637</td><td> 1.183974</td></tr>
	<tr><td>30.214595</td><td>31.149055</td><td>35.162583</td><td>34.973798</td><td>35.308077</td><td>34.814376</td><td>35.014091</td><td>26.699821</td><td>26.094639</td></tr>
	<tr><td>32.266133</td><td>24.629172</td><td>24.266737</td><td>22.514719</td><td>31.710218</td><td>35.212845</td><td>22.663506</td><td>18.210186</td><td>27.284313</td></tr>
	<tr><td>23.130484</td><td>21.403489</td><td>31.100670</td><td>31.830648</td><td>31.325406</td><td>31.622881</td><td>31.941189</td><td>17.946648</td><td>20.067772</td></tr>
	<tr><td>21.768145</td><td>32.941640</td><td>30.369523</td><td>28.210389</td><td>35.321865</td><td>27.567218</td><td>28.349570</td><td>29.263873</td><td>17.744246</td></tr>
	<tr><td>31.058407</td><td>32.076550</td><td>35.272635</td><td>34.427840</td><td>35.154832</td><td>34.292254</td><td>34.488305</td><td>25.946455</td><td>25.136552</td></tr>
	<tr><td>16.163062</td><td>22.309687</td><td>23.654234</td><td>19.539570</td><td>27.764840</td><td>19.810963</td><td>19.683332</td><td>12.142026</td><td> 8.683729</td></tr>
	<tr><td>26.331708</td><td>27.128958</td><td>32.232352</td><td>28.980514</td><td>31.793389</td><td>30.106267</td><td>29.107831</td><td>20.111331</td><td>19.866081</td></tr>
</tbody>
</table>
</dd>
	<dt>$design</dt>
		<dd><table class="table-responsive table-striped">
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
</dd>
</dl>




![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_47_2.png)


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
<font color ='#00bcd4'> In [23]: </font>

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
<font color ='#00bcd4'> In [24]: </font>

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
<font color ='#00bcd4'> In [25]: </font>

{% highlight R %}
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
{% endhighlight %}


2409


<br>
<font color ='#00bcd4'> In [26]: </font>

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
<font color ='#00bcd4'> In [27]: </font>

{% highlight R %}
vennDiagram(dt[,1:2], circle.col=c("mediumpurple", "midnightblue"))
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_55_0.png)


Venn diagram showing the number of genes DE in the comparison between basal
versus LP only (left), basal versus ML only (right), and the number of genes
that are DE in both comparisons (center). The number of genes that are not DE in
either comparison are marked in the bottom-right.

<br>
<font color ='#00bcd4'> In [28]: </font>

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
<font color ='#00bcd4'> In [29]: </font>

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
<font color ='#00bcd4'> In [37]: </font>

{% highlight R %}
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1],
       xlim=c(-8,13), col = c('mediumpurple', 'orange'))
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_61_0.png)


*Glimma* extends this functionality by providing an interactive mean-difference
plot via the `glMDPlot` function. The output of this function is an html page,
with summarised results in the left panel (similar to what is output by
`plotMD`), and the `log-CPM` values from individual samples in the right panel,
with a table of results below the plots. This interactive display allows the
user to search for particular genes based on their Gene symbol, which is not
possible in a static *R* plot. The `glMDPlot` function is not limited to mean-
difference plots, with a default version allowing a data frame to be passed with
the user able to select the columns of interest to plot in the left panel
(interactive plot  <a href="{{site.url}}{{site.baseurl}}/MDplot"> here </a>).

<br>
<font color ='#00bcd4'> In [31]: </font>

{% highlight R %}
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         id.column="ENTREZID", counts=x$counts, groups=group, launch=FALSE)
{% endhighlight %}

In the interactive plot the summary data (log-FCs
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
<font color ='#00bcd4'> In [32]: </font>

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


![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_65_1.png)


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
<font color ='#00bcd4'> In [33]: </font>

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
<font color ='#00bcd4'> In [34]: </font>

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
<font color ='#00bcd4'> In [35]: </font>

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
<font color ='#00bcd4'> In [36]: </font>

{% highlight R %}
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP,
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")
{% endhighlight %}


![png]({{ site.url}}{{ site.baseurl }}/notebooks/RNA-sequencing_files/RNA-sequencing_72_0.png)


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
