\name{GOSYNONYM}
\alias{GOSYNONYM}
\title{Map from GO synonyms to GO terms}
\description{
  GOSYNONYM is an R object that provides mapping
  from GO synonyms to GO terms
}

\details{
  [TODO: Put some details here]
}

\seealso{
  \itemize{
    \item \code{\link[AnnotationDbi]{AnnotationDb-class}} for use of
          the \code{select()} interface.
  }
}

\examples{
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.

## Bimap interface:
x <- GOSYNONYM
sample(x, 3)
# GO ID "GO:0009435" has a lot of synonyms
GOTERM[["GO:0009435"]]
# GO ID "GO:0006736" is a synonym of GO ID "GO:0009435"
GOID(GOSYNONYM[["GO:0006736"]])
}
\keyword{datasets}

