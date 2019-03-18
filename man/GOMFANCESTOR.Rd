\name{GOMFANCESTOR}
\alias{GOMFANCESTOR}
\title{Annotation of GO identifiers to their Molecular Function Ancestors} 
\description{
  This data set describes associations between GO molecular function (MF)
  terms and their ancestor MF terms, based on the directed acyclic
  graph (DAG) defined by the Gene Ontology Consortium. The format is an R
  object mapping the GO MF terms to all ancestor terms, where an
  ancestor term is a more general GO term that precedes
  the given GO term in the DAG (in other words, the parents, and all
  their parents, etc.).
}

\details{
  Each GO MF term is mapped to a vector of ancestor GO MF terms.

  Molecular function is defined as the tasks performed by individual
  gene products; examples are transcription factor and DNA helicase as
  defined by Gene Ontology Consortium.  
  
  Mappings were based on data provided: Gene Ontology  
  ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/  
  With a date stamp from the source of: 2018-Oct10

}

\references{
\url{http://www.geneontology.org/} and
\url{https://www.ncib.nlm.nih.gov/entrez/query.fcgi?db=gene} 
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
xx <- as.list(GOMFANCESTOR)
# Remove GO identifiers that do not have any ancestor
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # Get the ancestor GO identifiers for the first two elents of xx
  goids <- xx[1:2]
}
}
\keyword{datasets}

