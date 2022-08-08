# StaphCarrierTypes
Codes used to generate results for "Staphylococcus aureus Carrier Types are not Evidence of Population Heterogeneity"


The main markdown document is called getResults.rmd. If you download the repository and try knitting the .rmd, it should generate a document with figures (generating an HTML is fastest for me.) I have included Rdata files so you can (hopefully) generate all of the figures from the manuscript without running the optimization, MCMC, or simulation codes, which can take quite a while to run. The code chunks which do those things have been set to "eval=FALSE", so that they do not run if you knit the entire document. If you want to run simulations on your own, you can navigate to those specific code chunks.

The simulation codes towards the end of the document are set up to use parameter estimates from the earlier part of the document - so if you try running simulations and get an error, that might be why. 

The other most important files are the ones which define the log-likelihood models:

getL_ih_markdown.R

getL_2cr_markdown.R

getL_2s_markdown.R

getL_ucc_markdown.R

the one which loads in/formats the data:

dataLoader.R

and the data itself:

hcwDat05022019.csv

The codes to find Maximum Likelihood Estimates and generate Bayesian posterios are inside code chunks directly in the markdown document, although the code

MH_ensemble_markdown_gz.R 

is necessary to perform the affine-invariant ensemble sampling scheme described in the manuscript.

The codes

getSims.R

IHSim_gen.R 

2CRSim_gen.R 

UCCSim_gen.R 

2SSim_gen.R 

are needed to run simulations from the models.


The files in the respository which begin with

opt_

are the results I saved from optimizing the log-likelihood funtions, and these contain the Maximum Likelihood Estimates. At various places in the markdown file, I load in files like these to generate tables with the MLE under various assumptions. 

The files which begin with

jackknife_ 

are saved results from the jackknife scheme carried out in the markdown file.

The files which begin with

pmc_

are the output from the Monte Carlo sampling scheme which I used for the different models, and are loaded in to generate Bayesian posteriors.

The files which begin with

sims_

are output from simulations.
