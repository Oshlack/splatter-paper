Splatter paper
==============

Data and analysis for the paper "Splatter: Simulation of Single-cell RNA
sequencing data".

Directory structrure
--------------------

* `additional` - Supplementary figures and files
* `analysis` - Analysis files and output
* `data` - Input data files
* `figures` - Figures used in the main paper
* `output` - Additional intermediate files
* `R` - R functions used in analysis

Data
----

Data files used in the analysis are available in the `data.tar.gz` file. This
file should be extracted to `data` before attempting to run any of the analysis.

```{bash}
tar -xzvf data.tar.gz
```

After extraction the `data` directory will contain the following files:

* `datasets.txt` - Metadata about the various datasets
* `Camp.txt` - The Camp dataset
* `Engel.tsv` - The Engel dataset
* `Klein.csv` - The Klein dataset
* `Tung.txt` - The Tung dataset
* `Zeisel.txt` - The Zeisel dataset

Analysis
--------

The code for completing the analysis shown in the paper is provided as the
following Rmarkdown files in the analysis directory:

* `simulations.Rmd` - Some examples of Splat simulations.
* `datasets.Rmd` - Comparison of simulations based on various datasets.
* `clustering.Rmd` - Example evaluation of the SC3 clustering method.
* `supplementary.Rmd` - Supplementary figures.

There are two additional Rmd files which are rendered via `supplementary.Rmd`.
These are `additional_figures.Rmd` which combines the additional figures into a
single PDF and `sessionInfo.Rmd` which outputs the details of all the packages
used during the analysis.

Running the analysis files will produce figure files in the `figures` and
`additional` directories as well as data files in the `output` directory.

Please be aware that some of the analysis (particularly `datasets.Rmd`) requires
large amounts of resources (processing, memory, time) and may require slight
modifications to run in your environment.

R
---

This directory contains the following functions used in the analysis:

* `load_datasets.R`
  * `loadDataset` - Takes a row from `datasets.txt` and the path to the data
    files then returns an expression matrix for that dataset
* `simulate_datasets.R`
  * `simCompDataset` - Loads a dataset, estimates parameters for the various
    simulations, simulates new datasets and returns a comparison
* `utils.R`
  * `chrRound` - Rounds a number for presentation, converting it to a string
  * `logistic` - Implementation of the logistic function
  * `mcri.palettes` - List of colour palettes used for some plots
  * `mcriPalette` - Returns colour palette of particular size

