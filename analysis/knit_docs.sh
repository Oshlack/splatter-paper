#!/bin/bash
#PBS -N cellranger
#PBS -q batch
#PBS -l nodes=1:ppn=12
#PBS -l mem=64GB
#PBS -l walltime=20:00:00
#PBS -A bioinf
#PBS -m abe
#PBS -M luke.zappia@mcri.edu.au

#### Run:
module load R
module load pandoc
Rscript -e "library(rmarkdown)" \
    -e "rmarkdown::render('simulations.Rmd', 'html_document')"
    -e "rmarkdown::render('datasets.Rmd', 'html_document')"
    -e "rmarkdown::render('clustering.Rmd', 'html_document')"
    -e "rmarkdown::render('supplementary.Rmd', 'html_document')"
    -e "rmarkdown::render('sessionInfo.Rmd', 'html_document')"
