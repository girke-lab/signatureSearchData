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
qr <- query(ah, c("signatureSearchData", "taurefList"))
# taurefList <- ah[["AH69088"]]
}

\keyword{datasets}

