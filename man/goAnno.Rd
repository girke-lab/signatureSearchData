\name{goAnno}
\alias{goAnno}
\title{GO Term to Gene SYMBOLs Mapping Annotation} 
\description{
It is a \code{data.frame} containing the GO term ID to gene SYMBOLs mapping 
information used for \code{dtnetplot} function in \code{signatureSearch} package
to get the annotated gene labels set if the 'set' argument is a GO term ID.
}

\examples{
library(ExperimentHub)
eh <- ExperimentHub()
goAnno <- eh["EH3229"]
# goAnno <- eh[["EH3229"]]
}
\keyword{datasets}

