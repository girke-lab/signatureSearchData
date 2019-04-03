\name{ES_NULL}
\alias{ES_NULL}
\title{Null distribution of Enrichment Score from \code{gess_lincs} method} 
\description{
ES null distribution is generated with 1000 random queries searching against 
the \code{lincs} database for computing nominal 
P-values of Enrichment Score in \code{gess_lincs} result from 
\code{signatureSearch} package.
}
\details{
For more details, please refer to Subramanian et al., 2017, Cell, A Next 
Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles
}

\references{
  \url{https://www.sciencedirect.com/science/article/pii/S0092867417313090}
}

\seealso{
  \itemize{
    \item \code{signatureSearch}
    \item \code{gess_lincs}
  }
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "ES_NULL"))
ES_NULL <- qr[["AH69089"]]
}

\keyword{datasets}

