---
layout: page
permalink: /local_install/
title: Setting up R and Jupyter notebooks locally
category: page
---


# Getting all set up

Before being able to run the Jupyter notebooks locally you will need to install some dependencies. The basic instructions to do so are described below. Please make sure to install the components in the order detailed in this guide.


## 1. Installing R
You can find information regarding the latest R version on the [CRAN website](https://www.r-project.org).

* MAC OS-X:
Get the R binary from https://cran.r-project.org/bin/macosx/. Once downloaded install from the .pkg file.

* Windows:
Download the installer from (https://cran.r-project.org/bin/windows/base/). Once downloaded run the executable.

## 2. Installing Anaconda and the Notebooks
We *highly* recommend that you install the [Anaconda distribution](https://docs.continuum.io/anaconda/install) (or [Miniconda](https://conda.io/miniconda.html) alternatively).

You can download and install Anaconda on Windows, OSX and Linux. To ensure that it's up to date, run (in a terminal)

```
conda update conda
conda update jupyter
```
**Experienced users**

If you already have Python installed and prefer not to install Anaconda you can install the notebooks via pip:

```
pip install jupyter
```
You have now a notebook server installed on your computer. If you want to run a notebook server you need to open a terminal and run

```
jupyter notebook
```
Alternatively, if you have the Anaconda navigator you can open an instance from there.

For more information on running the notebooks server visit: https://jupyter.readthedocs.io/en/latest/running.html#running

## 3. Installing the R Kernel
This course will be using the [IRKernel](https://github.com/IRkernel/IRkernel). This has not been made available as a package from CRAN (yet). So in order to install the Kernel you need to install it via the `devtools` package (see [here](https://irkernel.github.io/installation/)):

```R
install.packages('devtools')
devtools::install_github('IRkernel/IRkernel')
# or devtools::install_local('IRkernel-master.tar.gz')
IRkernel::installspec()  # to register the kernel in the current R installation
```

Note you need to do this from an R console (**do not** use R studio if you have this installed).


## 4. Multiple packages
You will need to install various packages for this course. All of them are available on CRAN. Thus you can install them from a R console by typing

```R
install.packages('package_name')
```
Make sure you have all of the following packages:

* multtest (note that for the latest versions of R you need to install this from an R console, see the [Bioconductor website](https://bioconductor.org/packages/release/bioc/html/multtest.html)

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
```
* xlsx
* gdata
* ape


<br>
<a href="{{site.url}}{{site.baseurl}}/index.html" class="float">
<i class="fa fa-home my-float"></i>
</a>
