\name{lincs_expr}
\alias{lincs_expr}
\title{lincs expression signature database} 
\description{
The \code{lincs_expr} database contains mean gene expression values 
of 5,925 compound treatment samples in 30 cells (38,824 signatures in total). 
}
\details{
The original LINCS Level 3 expression database was downloaded at GEO 
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742}. 
The downloaded database in `gctx` format was processed and saved as an HDF5 file
}

\references{
  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742} 
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("signatureSearchData", "lincs_expr"))
# lincs_expr <- ah[["AH69093"]]
}

\keyword{datasets}

