\name{GOMFPARENTS}
\alias{GOMFPARENTS}
\title{Annotation of GO Identifiers to their Molecular Function Parents} 
\description{
  This data set describes associations between GO molecular function (MF)
  terms and their direct parent MF terms, based on the directed acyclic
  graph (DAG) defined by the Gene Ontology Consortium. The format is an R
  object mapping the GO MF terms to all direct parent terms, where a
  direct parent term is a more general GO term that immediately precedes
  the given GO term in the DAG.
}

\details{
  Each GO MF term is mapped to a named vector of GO MF terms. The name
  associated with the parent term will be either \emph{isa}, \emph{hasa} or \emph{partof},
  where \emph{isa} indicates that the child term is a more specific version
  of the parent, and \emph{hasa} and \emph{partof} indicate that the
  child term is a part of the parent. For example, a telomere is part of
  a chromosome.
  
  Molecular function is defined as the tasks performed by individual
  gene products; examples are transcription factor and DNA helicase as
  defined by the Gene Ontology Consortium.
  
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
xx <- as.list(GOMFPARENTS)
# Remove GO identifiers that do not have any parent
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
    # Get the parent GO identifiers for the first elents of xx
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

