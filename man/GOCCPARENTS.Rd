\name{GOCCPARENTS}
\alias{GOCCPARENTS}
\title{Annotation of GO Identifiers to their Cellular Component Parents} 
\description{
  This data set describes associations between GO molecular function (CC)
  terms and their direct parent CC terms, based on the directed acyclic
  graph (DAG) defined by the Gene Ontology Consortium. The format is an R
  object mapping the GO CC terms to all direct parent terms, where a
  direct parent term is a more general GO term that immediately precedes
  the given GO term in the DAG.
}
\details{
  Each GO CC term is mapped to a named vector of GO CC terms. The name
  associated with the parent term will be either \emph{isa}, \emph{hasa} or \emph{partof},
  where \emph{isa} indicates that the child term is a more specific version
  of the parent, and \emph{hasa} and \emph{partof} indicate that the
  child term is a part of the parent. For example, a telomere is part of
  a chromosome.
  
  Cellular component is defined as the subcellular structures,
  locations, and macromolecular complexes; examples include nucleus,
  telomere, and origin recognition complex as defined b y Gene Ontology
  Consortium.

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
xx <- as.list(GOCCPARENTS)
# Remove GO IDs that do not have any parent
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
   goids <- xx[[1]]
   # Find out the GO terms for the first parent go ID
   GOID(GOTERM[[goids[1]]])
   Term(GOTERM[[goids[1]]])
   Synonym(GOTERM[[goids[1]]])
   Secondary(GOTERM[[goids[1]]])
   Definition(GOTERM[[goids[1]]])
   Ontology(GOTERM[[goids[1]]])
}
}
\keyword{datasets}

