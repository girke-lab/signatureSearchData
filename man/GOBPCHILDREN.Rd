\name{GOBPCHILDREN}
\alias{GOBPCHILDREN}
\title{Annotation of GO Identifiers to their Biological Process Children} 
\description{
This data set describes associations between GO molecular function (BP)
  terms and their direct children BP terms, based on the directed acyclic
  graph (DAG) defined by the Gene Ontology Consortium. The format is an R
  object mapping the GO BP terms to all direct children terms, where a
  direct child term is a more specific GO term that is immediately preceded
  by the given GO term in the DAG.
}
\details{
   Each GO BP term is mapped to a vector of children GO BP terms.

  Biological process is defined as the broad biological goals, such as
  mitosis or purine metabolism, that are accomplished by ordered
  assemblies of molecular functions as defined by Gene Ontology
  Consortium.
  
  Mappings were based on data provided by: Gene Ontology  
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
xx <- as.list(GOBPCHILDREN)
# Remove GO IDs that do not have any children
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
    # Get the parent GO IDs for the first elents of xx
    goids <- xx[[1]]
    # Find out the GO terms for the first parent goid
    GOID(GOTERM[[goids[1]]])
    Term(GOTERM[[goids[1]]])
    Synonym(GOTERM[[goids[1]]])
    Secondary(GOTERM[[goids[1]]])
    Definition(GOTERM[[goids[1]]])
    Ontology(GOTERM[[goids[1]]])
}
}
\keyword{datasets}

