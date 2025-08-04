# genoPlotR

<p align="left"><img src="https://genoplotr.r-forge.r-project.org/img/genoplotr_logo.png" /></p>

genoPlotR is an R package made for generating reproducible, publication-grade 
graphics of gene and genome maps. It allows the user to read from various 
formats such as GenBank and BLAST results, as well as home-made tabular files.

<p align="left"><img src="https://genoplotr.r-forge.r-project.org/img/barto_multiseg.png" /></p>

# Installation

To install the latest stable version of genoPlotR:

    install.packages("genoPlotR")

To install the latest development version of genoPlotR from GitHub, the 
devtools package has to be installed first:

    install.packages("devtools")
    
If devtools is already installed, genoPlotR can be installed as follows:

    devtools::install_github("MVPuijk/genoPlotR")
    
Finally, if you want to run the command-line wrapper script `run_genoplot.R`, 
you must also install the optparse package:

    install.packages("optparse")
