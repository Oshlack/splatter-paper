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

Data
----

Data files used in the analysis are available in the `data.tar.gz` file. This
file should be extracted to `data` before attempting to run any of the analysis.

```{bash}
tar -xzvf data.tar.gz
```

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

Please be aware that some of the analysis (particularly `datasets.Rmd`) requires
large amounts of resources (processing, memory, time) and may require slight
modifications to run in your environment.

