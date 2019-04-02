\name{cmap_expr}
\alias{cmap_expr}
\title{cmap expression signature database}
\description{
For correlation based GESS methods in \code{\link{signatureSearch}} package, 
the signature database containing gene expression values from drug treatment 
samples could also be used to search for similarity. 

The `cmap_expr` database represents mean expression values of 1,309 drug
treatment samples in 4 cells (3,587 signatures in total).

It is stored as HDF5 backed `SummarizedExperiment` object.
}
\details{
The same as the \code{\link{cmap}} database, the resource data was downloaded 
from the Connectivity Map project 
\url{https://portals.broadinstitute.org/cmap/}.
The expression intensity data was rma normalized. Then, the mean expression 
values of drug treatment samples at different concentration, duration was used
as expression values of that drug in a specific cell. The expression matrix
was saved as HDF5-based SummarizedExperiment object (cmap_expr database).
}

\seealso{
  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "cmap"))
qr
# se <- load_sigdb("AH69077","AH69078","cmap_expr")
}

\keyword{datasets}

