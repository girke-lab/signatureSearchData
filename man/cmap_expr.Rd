\name{cmap_expr}
\alias{cmap_expr}
\title{cmap expression signature database}
\description{
For correlation based GESS methods in the \code{signatureSearch} software package, 
the signature database containing gene expression values of drug treatment 
samples could also be used to search for similarity. 

The `cmap_expr` database represents mean expression values of 1,309 drug
treatment samples in 4 cells (3,587 signatures in total).
It is stored in an HDF5 file.
}
\details{
The same as the \code{\link{cmap}} database, the resource data was downloaded 
from the Connectivity Map project 
\url{https://portals.broadinstitute.org/cmap/}.
The expression intensity data was normalized by MAS5 method. Then, the mean expression 
values of drug treatment samples at different concentration, duration was used
as expression values of that drug treatment in a cell. The expression matrix
was written in an HDF5 file (cmap_expr database).
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("signatureSearchData", "cmap_expr"))
# cmap_expr <- ah[["AH69091"]]
}

\keyword{datasets}

