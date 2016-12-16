As scRNA-seq data has become available there has been a rapid development of new bioinformatics tools
attempting to unlock it's potential. Currently there are at least 80 available software packages that
have been designed specifically for the analysis of scRNA-seq data, the majority of which have been published
in peer-reviewed journals or as preprints. A table of scRNA-seq software is available at [https://goo.gl/4wcVwn](https://goo.gl/4wcVwn). The focus of these packages is often different from those
designed for the analysis of a bulk RNA-seq experiment. In a bulk experiment the groups of samples are known
and the task is often to test for genes that are differentially expressed (DE) between two or more groups. In contrast the groups
in a single-cell experiment are usually unknown and the analysis is often more exploratory. Many of the
current packages focus on the task of assigning cells to groups before applying more traditional differential expression testing. This approach is taken by tools such as SC3 \cite{Kiselev2016-fa}, CIDR \cite{Lin2016-yu} and Seurat \cite{Satija_2015} and makes sense for a sample with a defined set of mature cell types. In the developmental setting, for example, where stem cells are differentiating into mature cells, it may be more appropriate
to order cells along a continuous trajectory from one cell type to another. Tools such as Monocle \cite{Trapnell_2014}, CellTree \cite{duVerle_2016} and 
Sincell \cite{Julia2015-zc} take this approach, ordering cells along a path then looking for patterns in the changes of gene expression.