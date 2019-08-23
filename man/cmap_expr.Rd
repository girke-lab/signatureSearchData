\name{cmap_expr}
\alias{cmap_expr}
\title{CMap Intensity Signature Database}
\description{
To search CMap2 with signatureSearch's correlation based GESS methods (gess_cor),
normalized gene expression values (here intensities) are required where the 
biological replicate information has been collapsed to mean values. For this, 
the \code{cmap_expr} database has been created from CEL files, which are the 
raw data of the Affymetrix technology. To obtain normalized expression data,
the CEL files were downloaded from the CMap project site, and then processed 
with the MAS5 algorithm. Gene level expression data was generated the same way 
as described in help file of \code{cmap_rank}. Next, the gene expression values 
for different concentrations and treatment times of each compound and cell 
were averaged. Subsequently, the expression matrix was saved to an HDF5 file, 
referred to as the cmap_expr database. It represents mean expression values of 
12,403 genes for 1,309 compound treatments in up to 5 cells (3,587 signatures in total).
}
\details{
The cmap_expr database can be downloaded as HDF5 file from Bioconductor’s 
AnnotationHub as shown in the Example section. 

The loaded \code{cmap_expr} data object is generated from the raw CEL files
downloaded from the CMap project site.
For documentation and code of generating the cmap_expr databases from sources, 
please refer to the vignette of this package by running 
\code{browseVignettes("signatureSearchData")} in R.
}

\references{
CMap project site: https://portals.broadinstitute.org/cmap
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("signatureSearchData", "cmap_expr"))
# cmap_expr_path <- ah[["AH69091"]]
}

\keyword{datasets}

