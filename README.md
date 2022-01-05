# _signatureSearchData_: Reference Data for Gene Expression Signature Searching

# Introduction

The `signatureSearchData` package provides access to the reference data used by
the associated `signatureSearch` software package. The latter allows to search 
with a query gene expression signature (GES) a database of treatment GESs to 
identify cellular states sharing similar expression responses (connections). This 
way one can identify drugs or gene knockouts that induce expression phenotypes 
similar to a sample of interest. The resulting associations may lead to novel 
functional insights how perturbagens of interest interact with biological systems. 

Currently, `signatureSearchData` includes GES data from the CMap (Connectivity
Map) and LINCS (Library of Network-Based Cellular Signatures) projects that are
largely based on drug and genetic perturbation experiments performed on
variable numbers of human cell lines [@Lamb2006-du; @Subramanian2017-fu]. In
`signatureSearchData` these data sets have been preprocessed to be compatible
with the different gene expression signature search (GESS) algorithms
implemented in `signatureSearch`. The preprocessed data types include but are
not limited to normalized gene expression values (_e.g._ intensity values), log
fold changes (LFC) and Z-scores, p-values or FDRs of differentially expressed
genes (DEGs), rankings based on selected preprocessing routines or sets of top
up/down-regulated DEGs. 

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
and treatment times to ~70 human non-cancer (normal) and cancer cell lines.
Additionally, it includes several thousand genetic perturbagens composed of
gene knockdown and over-expression experiments. 

In 2020, the LINCS 2017 database is expanded to the beta release,
here refer to as LINCS2. It contains >80k perturbations and >200 cell lines and 
over 3M gene expression profiles. This represents roughly a 3-fold expansion on 
the LINCS 2017 database and notable new subsets of data include CRISPSR knockout 
of >5k genes and hematopoietic and non-cancer cell models. 
The datasets can be accessed at https://clue.io/releases/data-dashboard.

The data structures and search algorithms used by `signatureSearch` and
`signatureSearchData` are designed to work with most genome-wide expression
data including hybridization-based methods, such as Affymetrix or L1000, as
well as sequencing-based methods, such as RNA-Seq. Currently,
`signatureSearchData` does not include preconfigured RNA-Seq reference data mainly 
due to the lack of large-scale perturbation studies (_e.g._ drug-based) available 
in the public domain that are based on RNA-Seq. This situation may change in 
the near future once the technology has become more affordable for this purpose. 

# Vignette
The vignette of this package is available at [here](https://www.bioconductor.org/packages/release/data/experiment/vignettes/signatureSearchData/inst/doc/signatureSearchData.html)

# Install and Load Package

`signatureSearchData` is a R/Bioconductor package and can be installed using 
`BiocManager::install()`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("signatureSearchData", version = "3.9")
```

To obtain the most recent updates immediately, one can install it directly from GitHub as follows.
```r
devtools::install_github("yduan004/signatureSearchData")
```
After the package is installed, it can be loaded into an R session as follows.
```r
library(signatureSearchData)
```
For a detailed description of loading and/or generating signature databases,
please refer to the vignette of this package by running
`browseVignettes("signatureSearchData")` in an R session.
