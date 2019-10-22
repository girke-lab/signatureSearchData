\name{GO_DATA}
\alias{GO_DATA}
\title{GO Annotation Environment} 
\description{
It is an RDS file storing the GO annotation environment (GO term ID to 
gene SYMBOLs mapping, gene SYMBOL to annotated GO terms, GO term ID to 
its description, GO term ID to its ontology) 
used for \code{tsea_*} functions in \code{signatureSearch} package with GO 
enrichment analysis. This environment is saved to an RDS file to accelerate the 
analysis speed by avoiding building this environment from scratch every time 
running the \code{tsea} functions.
}

\seealso{
    \code{tsea}
}

\examples{
library(ExperimentHub)
eh <- ExperimentHub()
qr <- query(eh, c("signatureSearchData", "GO_DATA"))
go_data <- eh["EH3231"]
# go_data <- eh[["EH3231"]]
}

\keyword{datasets}

