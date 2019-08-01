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
please refer to the vignette of the \code{signatureSearch} package by
running \code{browseVignettes("signatureSearch")} in your R session.
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

