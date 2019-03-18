\name{GOOBSOLETE}
\alias{GOOBSOLETE}
\title{Annotation of GO identifiers by terms defined by Gene Ontology
  Consortium and their status are obsolete}
\description{
This is an R object mapping GO identifiers to the specific
terms in defined by Gene Ontology Consortium and their definition are obsolete
}
\details{
  All the obsolete GO terms that are collected in this index will no longer exist 
  in other mapping objects.

  Mappings were based on data provided by: Gene Ontology  
  ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/  
  With a date stamp from the source of: 2018-Oct10

}
\references{
  \url{https://www.ncib.nlm.nih.gov/LocusLink} 
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
# Convert the object to a list
xx <- as.list(GOTERM)
if(length(xx) > 0){
    # Get the TERMS for the first elent of xx
    GOID(xx[[1]])
    Ontology(xx[[1]])
}
}
\keyword{datasets}

