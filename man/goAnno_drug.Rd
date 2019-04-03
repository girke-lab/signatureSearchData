\name{goAnno_drug}
\alias{goAnno_drug}
\title{Annotation of GO ID to drug name mappings} 
\description{
It is an intermediate file storing the GO ID to drug name mapping information
for \code{dsea_*} functions in \code{signatureSearch} package.
}

\details{
The GO ID-GENE SYMBOL mappings were converted to GO ID-drug mappings via
drug-target links in DrugBank, CLUE and STITCH databases
}
\seealso{
  \itemize{
    \item \code{signatureSearch}
    \item \code{dsea_hyperG}
    \item \code{dsea_GSEA}
  }
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "goAnno"))
# goAnno_drug <- qr[["AH69085"]]
}

\keyword{datasets}

