#!/bin/bash
#PBS -N knit_docs
#PBS -q batch
#PBS -l nodes=1:ppn=12
#PBS -l mem=200GB
#PBS -l walltime=100:00:00
#PBS -A bioinf
#PBS -m abe
#PBS -M luke.zappia@mcri.edu.au

#### Run:
module load R
module load pandoc

cd /group/bioi1/luke/analysis/splatter-paper/analysis

Rscript -e "library(rmarkdown)" \
    -e "rmarkdown::render('datasets.Rmd', 'html_document')"
