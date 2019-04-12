\name{goAnno}
\alias{goAnno}
\title{Annotation of GO ID to gene SYMBOL mappings} 
\description{
 It is an intermediate file storing the GO ID to gene SYMBOL mapping 
 information for \code{dtnetplot} function in \code{signatureSearch} package
}

\seealso{
    \code{dtnetplot}
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "goAnno"))
goAnno <- qr[["AH69084"]]
}
\keyword{datasets}

