\name{taurefList}
\alias{taurefList}
\title{Reference list for computing tau score of \code{gess_lincs} method} 
\description{
It uses all signatures in the LINCS reference database to query against 
itself as Qref to compute tau score of \code{gess_lincs} 
method in \code{signatureSearch} package. Tau score compares observed enrichment
score to all others in Qref. It represents the percentage of reference queries 
with a lower |NCS| than |NCSq,r|, adjusted to retain the sign of NCSq,r. 
NCSq,r is the normalized connectivity score for signature r relative to 
query q. A tau of 90 indicates that only 10 percent of reference perturbations 
showed stronger connectivity to the query. 
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
qr <- query(ah, c("signatureSearchData", "taurefList"))
# taurefList <- qr[["AH69088"]]
}

\keyword{datasets}

