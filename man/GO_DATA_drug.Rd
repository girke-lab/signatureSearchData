\name{GO_DATA_drug}
\alias{GO_DATA_drug}
\title{GO to Drug Annotation Environment} 
\description{
It is an RDS file storing the GO annotation environment where the gene SYMBOLs 
are mapped to drug names via drug-target interactions. The GO to drug annotation
environment contains GO term ID to drug names, drug name to annotated 
GO terms, GO term ID to its description and GO term ID to its ontology mappings.
This environment is used for \code{dsea_*} methods in the \code{signatureSearch} 
package to do GO enrichment analysis directly with test drug set.
This environment is stored in an RDS file to accelerate the speed by avoiding 
building this environment from scratch every time running the \code{dsea} functions. 
}

\seealso{
    \code{dsea}
}

\examples{
library(ExperimentHub)
eh <- ExperimentHub()
qr <- query(eh, c("signatureSearchData", "GO_DATA"))
go_data_drug <- eh["EH3232"]
# go_data_drug <- eh[["EH3232"]]
}

\keyword{datasets}

