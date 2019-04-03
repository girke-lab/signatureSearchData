\name{goAnno}
\alias{goAnno}
\title{Annotation of GO ID to GENE SYMBOL mappings} 
\description{
 It is an intermediate file storing the GO ID to GENE SYMBOLE mapping 
 information for \code{dtnetplot} function in \code{signatureSearch} package
}

\seealso{
  \itemize{
    \item \code{signatureSearch}
    \item \code{dtnetplot}
  }
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "goAnno"))
goAnno <- qr[["AH69084"]]
}
\keyword{datasets}

