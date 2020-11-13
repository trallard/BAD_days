# Bioinformatics Awareness Days (BAD).

Originally developed by [Luisa Cutillo](https://www.sheffield.ac.uk/neuroscience/staff/cutillo).

The Bioinformatics Awareness Days are days devoted to Bioinformatics. As a Marie S. Curie postdoctoral fellow @ SITRaN and Computer Science Departments at the University of Sheffield, I organized such an event to give my practical contribution in Bioinformatics, although I have to state that all this started in Naples at the Telethon Institute of Genetics and Medicine, with the great Bionformatics Core.

Together with [Tania Allard](https://github.com/trallard) and [Mike Croucher](https:/github.com/mikecroucher) at the DCS, we decided to publicly disseminate this material to all the interested scientific community. The sessions are self contained and a full run should last at most 2 hours.

**All of the sessions material is now also contained in a website based around the Jupyter notebooks: [BAD Days website](http://bitsandchips.me/BAD_days/).**

## What is R?

R is a free, open-source programming language that has very strong support for statistics. It was originally developed as an open source implementation of the [S Programming language](https://en.wikipedia.org/wiki/S_(programming_language)).
It is used extensively in research and industry for areas such as data analysis, statistics, machine learning, bioinformatics, simulation, linguistics and much more.

With over [8000 freely available add-on packages](https://cran.r-project.org/web/packages/) that provide extensive additional functionality, R will probably have something that can help your research.

Don't just take our word for it though -- here's what others have to say

* [Why use R? Five reasons](http://www.econometricsbysimulation.com/2014/03/why-use-r-five-reasons.html) - From the 'Econometrics By Simulation' blog.

# About the course

The core of this tutorial is built around the [Jupyter Notebooks](https://jupyter.org), an interactive computational environment run in a web browser.

In this particular case we will be using the R programming language for the data exploration/analysis.

---

# Working with the material
The simplest way to enjoy these materials is to view each lesson online [here](http://bitsandchips.me/BAD_days/). In such a case, you could follow along and write your own version of the programs to run with your local R installation.

You can also **download/fork** the whole collection of Notebooks from the GitHub repository or using the download button within the BAD days website. That way you can use your local copy of the materials and follow along the course and make your own modifications/extensions.


# Getting all set up

Before continuing with the material you will need to install some dependencies. The basic instructions to do so are described below. Please make sure to install the components in the order detailed in this guide.

## Working locally

### 1. Installing R
You can find information regarding the latest R version on the [CRAN website](https://www.r-project.org).

* MAC OS-X:
Get the R binary from https://cran.r-project.org/bin/macosx/. Once downloaded install from the .pkg file.

* Windows:
Download the installer from https://cran.r-project.org/bin/windows/base/. Once downloaded run the executable.

### 2. Installing Anaconda and the Notebooks
We *highly* recommend that you install the [Anaconda distribution](https://docs.continuum.io/anaconda/install) (or [Miniconda](https://conda.io/miniconda.html) alternatively).

You can download and install Anaconda on Windows, OSX and Linux. To ensure that it's up to date, run (in a terminal)
```
conda update conda
conda update jupyter
```
**Experienced users**
If you already have Python installed and prefer not to install Anaconda you can install the notebooks via pip:
```bash
pip install jupyter
```
You have now a notebook server installed on your computer. If you want to run a notebook server you need to open a terminal and run
```
jupyter notebook
```
Alternatively, if you have the Anaconda navigator you can open an instance from there.

For more information on running the notebooks server visit: https://jupyter.readthedocs.io/en/latest/running.html#running

### 3. Installing the R Kernel
This course will be using the [IRKernel](https://github.com/IRkernel/IRkernel). This has not been made available as a package from CRAN (yet). So in order to install the Kernel you need to install it via the `devtools` package (see [here](https://irkernel.github.io/installation/)):
```R
install.packages('devtools')
devtools::install_github('IRkernel/IRkernel')
# or devtools::install_local('IRkernel-master.tar.gz')
IRkernel::installspec()  # to register the kernel in the current R installation
```

Note you need to do this from an R console (**do not** use R studio if you have this installed).


### 4. Multiple packages
You will need to install various packages for this course. All of them are available on CRAN. Thus you can install them from a R console by typing
```R
install.packages('package_name')
```
Make sure you have all of the following packages:
* multtest (note that for the latest versions of R you need to install this from an R console, see the [Bioconductor website](https://bioconductor.org/packages/release/bioc/html/multtest.html))
```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
```
* xlsx
* gdata
* ape
* ggplot2
* R.utils
* Mus.musculus
* Glimma
* limma


## Working on the cloud
Some notebooks have been added to a library within [Microsoft Azure Notebooks](https://notebooks.azure.com/?WT.mc_id=academic-0000-taallard). Whenever you see the cloud action button (see figure above) you can launch the corresponding notebook in Microsoft Azure for free (see the website).


Some of the advantages of this method is that you do not need to install any packages/dependencies in your local computer.


Follow the first steps to work with the notebooks in the cloud:

1. Create a [Microsoft Azure Notebooks](https://notebooks.azure.com/?WT.mc_id=academic-0000-taallard). Head to the website and click on sign-in (upper right corner). Use your university email address and set a password.
2. Make sure you are logged into your account.
3. Go to the Notebook in this website that you want to run in the cloud and click on the Cloud icon located at the bottom right corner of your screen.
4. This should take you to the library where the notebooks are saved to.
5. Click on the notebook's line (it should highlight the whole line)![azure2](./images/azure2.png).
6. Click on the clone icon and create a new library to which this will be cloned. This will allow you to have your own copy of the notebook.
