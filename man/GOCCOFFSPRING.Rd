\name{GOCCOFFSPRING}
\alias{GOCCOFFSPRING}
\title{Annotation of GO Identifiers to their Cellular Component Offspring} 
\description{
  This data set describes associations between GO molecular function (MF)
  terms and their offspring MF terms, based on the directed acyclic
  graph (DAG) defined by the Gene Ontology Consortium. The format is an R
  object mapping the GO MF terms to all offspring terms, where an
  offspring term is a more specific GO term that is preceded
  by the given GO term in the DAG (in other words, the children and all
  their children, etc.).
}
\details{  
  Each GO CC term is mapped to a vector of offspring GO MF terms.
  
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
xx <- as.list(GOCCOFFSPRING)
# Remove GO identifiers that do not have any offspring
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # Get the offspring GO identifiers for the first two elents of xx
  goidentifiers <- xx[1:2]
}
}
\keyword{datasets}

