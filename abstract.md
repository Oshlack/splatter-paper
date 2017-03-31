**Motivation**
As single-cell RNA-seq (scRNA-seq) technology has rapidly developed so have methods for analysing it. Many of these methods have been tested and developed using simulated datasets. While this is a valid and useful approach there are several ways in which it can be improved. Many of the current simulations are not well documented, code may not be available and they don't demonstrate similarity to real data.

**Results**
We present the Splatter package for simple simulation of single-cell RNA-seq data. Splatter is a Bioconductor R package that aims to provide a consistent interface for multiple scRNA-seq simulation methods that are well-documented and easy to use. Our own simulation, Splat, is based on a gamma-poisson distribution incorporating high-expression outlier genes, defined library sizes, a mean-variance trend and technical dropout, allowing users to simulate single populations of cells, populations with multiple cell types or differentiation paths. Splatter also makes it easy to compare real and simulated datasets and we provide a short example of how the Splat simulation can be used to evaluate an analysis method.

**Availability and Implementation**
The Splatter R package, along with installation instructions and a vignette describing it's use, is available through Bioconductor at [https://www.bioconductor.org/packages/splatter/](https://www.bioconductor.org/packages/splatter/). The code is being developed on Github at [https://github.com/Oshlack/splatter](https://github.com/Oshlack/splatter).

**Contact**
_Full email address to be given, preferably an institution email address._

**Supplementary information**
_Links to additional figures/data available on a web site, pr reference to online-only supplementary data available at the journal's web site._