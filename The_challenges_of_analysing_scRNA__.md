The challenges of analysing scRNA-seq are also different to those of a bulk RNA-seq experiment. As there is such a small
amount of starting material in an individual cell, and we sequence many more samples, the sequencing
depth is lower in an scRNA-seq experiment. Related to this is the large number of genes for which
there is no measured expression. Some of this is due to the biology we wish to study. For example we expect
different cell types to express different genes, but there are additional factors such as stage in the
cell cycle, transcriptional bursting and environmental interactions which cause differences in expression
between cells performing the same function. There are also technical effects to consider. Due to the low
sequencing depth it is difficult to accurately measure expression in lowly expressed genes. Furthermore existing
protocols may not reliably capture all the RNA present, resulting in "dropout" events where a gene is expressed in
a sample but not observed in the sequencing data. Additionally it is impossible to replicate an
individual cell and it may not be possible to capture all the cell types of interest. Quality control
of damaged, missing, multiple or otherwise low-quality cells is also an important consideration.