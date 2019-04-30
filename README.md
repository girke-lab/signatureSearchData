# _signatureSearchData_ - Reference Data for Large-Scale Gene Expression Signature Searching

# Introduction

The `signatureSearchData` package provides access to the reference data used by
the associated `signatureSearch` software package. The latter allows to search 
with a query gene expression signature a database of treatment signatures to 
identify cellular states with similar expression responses (connections). This 
way one can identify drugs or gene knockouts that induce expression phenotypes 
similar to a sample of interest. The resulting associations may lead to novel 
functional insights how perturbagens of interest interact with biological systems. 

Currently, `signatureSearchData` includes gene expression data from the CMap 
(Connectivity Map) and LINCS (Library of Network-Based Cellular Signatures) 
projects that are largely based on drug and
genetic perturbation experiments performed on variable numbers of human cell
lines. In `signatureSearchData` these data sets 
are preprocessed as needed by the gene expression signature search (GESS)
algorithms implemented and used by `signatureSearch`. The preprocessed data
types include but are not limited to normalized gene expression values 
(_e.g._ intensity values), log fold changes (LFC) and Z-scores, P or FDR values 
for DEG calls, rankings based on selected preprocessing routines or sets of 
top up/down-regulated DEGs. 

The CMap data were downloaded from the [CMap project
site](https://portals.broadinstitute.org/cmap/) (Version build02). The latter is
a collection of over 7,000 gene expression profiles (signatures) obtained
from perturbation experiments with 1,309 drug-like small molecules on five
human cancer cell lines. The Affymetrix Gene Chip technology was used to
generate the CMAP2 data set. 

In 2017, the LINCS Consortium generated a similar but much larger data set where
the total number of gene expression signatures was scaled up to over one
million. This was achieved by switching to a much more cost effective gene
expression profiling technology called L1000 assay [@Peck2006-rf;
@Edgar2002-di]. The current set of perturbations covered by the LINCS data set
includes 19,811 drug-like small molecules applied at variable concentrations
and treatment times to 70 human non-cancer (normal) and cancer cell lines.
Additionally, it includes several thousand genetic perturbagens composed of
gene knockdown and over-expression experiments. 

The data structures and search algorithms used by `signatureSearch` and
`signatureSearchData` are designed to work with most genome-wide expression
data including hybridization-based methods, such as Affymetrix or L1000, as
well as sequencing-based methods, such as RNA-Seq. Currently,
`signatureSearchData` does not include RNA-Seq data mainly due to the lack
of large-scale perturbation studies (_e.g._ drug-based) available in the public
domain that are based on RNA-Seq. This situation may change in the near future
once the cost of RNA-Seq becomes more cost effective for this purpose. 

# Install and load `signatureSearchData`

`signatureSearch` is a R/Bioconductor package and can be installed using 
`BiocManager::install()`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("signatureSearchData", version = "3.9")
```
Or it can be directly installed from GitHub by
```r
devtools::install_github("yduan004/signatureSearchData")
```
After the package is installed, it can be loaded into an R session as follows.
```r
library(signatureSearchData)
```
For detailed description of generating the CMAP/LINCS signature databases, please refer to the vignette 
of this package by running `browseVignettes("signatureSearchData")` in R session.
