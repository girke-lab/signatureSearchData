\name{ES_NULL}
\alias{ES_NULL}
\title{Null distribution of WTCS score of \code{gess_lincs} method} 
\description{
WTCS null distribution is generated with 1000 random queries searching against 
the \code{lincs} database for computing nominal 
p-values of WTCS score in the \code{gess_lincs} result from 
\code{signatureSearch} package.
}

\references{
For detailed description of the \code{gess_lincs} method and scores, 
please refer to: Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D., 
Natoli, T. E., Lu, X., … Golub, T. R. (2017). A Next Generation 
Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell, 
171(6), 1437–1452.e17. \url{https://doi.org/10.1016/j.cell.2017.10.049}
}

\seealso{
   \code{gess_lincs}
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "ES_NULL"))
ES_NULL <- ah[["AH69089"]]
}

\keyword{datasets}

