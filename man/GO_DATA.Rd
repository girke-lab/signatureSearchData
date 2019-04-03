\name{GO_DATA}
\alias{GO_DATA}
\title{GO annotation environment} 
\description{
It is an intermediate file storing the GO annotation environment (GO term to 
gene sets, gene to GO terms, GO term ID to description, GO ID to ontology) 
used for \code{tsea_*} methods in \code{signatureSearch} package to do GO 
enrichment analysis to accelarate the speed by avoiding building this 
environment every time running the tsea functions 
}

\seealso{
  \itemize{
    \item \code{signatureSearch}
    \item \code{tsea}
  }
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "GO_DATA"))
#GO_DATA <- qr[["AH69086"]]
}

\keyword{datasets}

