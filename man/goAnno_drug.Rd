\name{goAnno_drug}
\alias{goAnno_drug}
\title{GO Term ID to Drug Name Mapping Annotation} 
\description{
It is a \code{data.frame} containing the GO term ID to drug name mapping 
information used for \code{dsea_*} functions in \code{signatureSearch} package.
}

\details{
The original GO annotation system contains GO term ID to gene SYMBOL id mappings,
the latter are converted to drug names mappings via
drug to target protein-encoding genes links in DrugBank, CLUE and STITCH databases.
}
\seealso{
  \itemize{
    \item \code{dsea_hyperG}
    \item \code{dsea_GSEA}
  }
}

\examples{
library(ExperimentHub)
eh <- ExperimentHub()
goAnno_drug <- eh["EH3230"]
# goAnno_drug <- eh[["EH3230"]]
}

\keyword{datasets}

