\name{signatureSearch_data}
\alias{signatureSearch_data}
\title{Annotation data package for signatureSearch}
\description{
This is an annotation data package for \code{signatureSearch} softeware packge. 
It stores the CMAP/LINCS reference databases used for GESS methods in the
\code{signatureSearch} software package as well as other intermediate datasets used for
functionalities in the \code{signatureSearch} package, 
such as getting targets of query drugs in DrugBank, CLUE and STITCH databases, 
conducting GO enrichment analysis or computing p-values or Tau scores 
in the results of `gess_lincs` method.
}
\details{
The datasets in this package are:
  \itemize{
    \item \code{\link{cmap}} 
    \item \code{\link{cmap_expr}} 
    \item \code{\link{lincs}} 
    \item \code{\link{lincs_expr}} 
    \item \code{\link{dtlink_db_clue_sti}} 
    \item \code{\link{goAnno}} 
    \item \code{\link{goAnno_drug}} 
    \item \code{\link{GO_DATA}} 
    \item \code{\link{GO_DATA_drug}} 
    \item \code{\link{taurefList}} 
    \item \code{\link{ES_NULL}} 
  }
For documentations of generating the CMAP/LINCS databases, please refer to the 
vignette of this package by running `browseVignettes("signatureSearchData")` in R.
}

\seealso{
   \code{signatureSearch}
}

\examples{
library(AnnotationHub)
# Load cmap database
# se <- load_sigdb("AH69075","AH69076","cmap")
# Load dtlink SQLite database
#conn <- load_sqlite("AH69083")
# Load RDS files
ah <- AnnotationHub()
#goAnno <- ah[["AH69084"]]
}

\keyword{package}

