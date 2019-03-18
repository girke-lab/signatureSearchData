\name{GOTERM}
\alias{GOTERM}
\title{Annotation of GO Identifiers to GO Terms}
\description{
  This data set gives mappings between GO identifiers and their respective terms.
}
\details{
  Each GO identifier is mapped to a \code{GOTerms} object that has 6 slots:
   GOID: GO Identifier
   Term: The term for that GO id
   Synonym: Synonymous terms
   Secondary: Secondary terms that have been merged into this term
   Definition: Further definition of the GO term
   Ontology: One of MF - molecular function, BP - biological
    process, or CC - cellular component
 
  All the obsolete GO terms are under the nodes "obsolete molecular
  function" (GO:0008369), "obsolete cellular component" (GO id
  GO:0008370), and "obsolete biological process" (GO:0008371). Each of
  these GO identifiers has a group of GO identifiers as their direct children with GO
  terms that were defined by GO but are
  deprecated in the current build. These deprecated GO terms were
  appended by "(obsolete)" when the data package was built. 

  Mappings were based on data provided by: Gene Ontology  
  ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/  
  With a date stamp from the source of: 2018-Oct10

}
\references{
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
xx <- as.list(GOTERM)
if(length(xx) > 0){
    # Get the TERMS for the first elent of xx
    GOID(xx[[1]])
    Term(xx[[1]])
    Synonym(xx[[1]])
    Secondary(xx[[1]])
    Definition(xx[[1]])
    Ontology(xx[[1]])
}
}
\keyword{datasets}

