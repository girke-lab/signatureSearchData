\name{GOMAPCOUNTS} 

\alias{GOMAPCOUNTS}


\title{Number of mapped keys for the maps in package GO.db}

\description{
  GOMAPCOUNTS provides the "map count" (i.e. the
  count of mapped keys) for each map in package GO.db.
}

\details{  
  This "map count" information is precalculated and stored in the
  package annotation DB. This allows some quality control and is used by
  the \code{\link[AnnotationDbi:AnnDbPkg-checker]{checkMAPCOUNTS}}
  function defined in AnnotationDbi to compare and validate different
  methods (like \code{count.mappedkeys(x)} or
  \code{sum(!is.na(as.list(x)))}) for getting the "map count" of a given
  map.
}

\seealso{
  \itemize{
    \item \code{\link[AnnotationDbi:Bimap-keys]{mappedkeys}},
    \item \code{\link[AnnotationDbi:Bimap-keys]{count.mappedkeys}},
    \item \code{\link[AnnotationDbi:AnnDbPkg-checker]{checkMAPCOUNTS}}
    \item \code{\link[AnnotationDbi]{AnnotationDb-class}} for use of
          the \code{select()} interface.
  }
}

\examples{
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.

## Bimap interface:
GOMAPCOUNTS
mapnames <- names(GOMAPCOUNTS)
GOMAPCOUNTS[mapnames[1]]
x <- get(mapnames[1])
sum(!is.na(as.list(x)))
count.mappedkeys(x)   # much faster!

## Check the "map count" of all the maps in package GO.db
checkMAPCOUNTS("GO.db")
}

\keyword{datasets}

