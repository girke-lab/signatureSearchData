\name{lincs}
\alias{lincs}
\title{lincs signature database} 
\description{
The \code{lincs} database represents moderated z-scores from DE analysis of 
12,328 genes from 8,140 compound treatments in 30 cells 
(45,956 signatures in total). 
}
\details{
The original LINCS Level 5 moderated z-score database was downloaded at GEO
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742}. 
The downloaded database in `gctx` format was processed and saved as an HDF5 file
}

\references{
  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742} 
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("signatureSearchData", "lincs"))
# lincs <- ah[["AH69092"]]
}

\keyword{datasets}

