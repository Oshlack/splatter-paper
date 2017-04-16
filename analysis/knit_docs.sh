#!/bin/bash
#PBS -N knit_docs
#PBS -q batch
#PBS -l nodes=1:ppn=12
#PBS -l mem=200GB
#PBS -l walltime=20:00:00
#PBS -A bioinf
#PBS -m abe
#PBS -M luke.zappia@mcri.edu.au

#### Run:
module load R
module load pandoc

cd /group/bioi1/luke/analysis/Simplifying-simulation-of-single-cell-RNA-sequencing-with-Splatter/analysis

Rscript -e "library(rmarkdown)" \
    -e "rmarkdown::render('simulations.Rmd', 'html_document')" \
    -e "rmarkdown::render('datasets.Rmd', 'html_document')" \
    -e "rmarkdown::render('clustering.Rmd', 'html_document')" \
    -e "rmarkdown::render('supplementary.Rmd', 'html_document')" \
    -e "rmarkdown::render('sessionInfo.Rmd', 'html_document')"
