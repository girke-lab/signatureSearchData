\name{signatureSearchData}
\alias{signatureSearchData}
\title{Annotation data package for signatureSearch}
\description{
This is an annotation data package for \code{signatureSearch} software package. 
It stores the CMAP/LINCS reference databases used for GESS methods in the
\code{signatureSearch} package as well as other intermediate datasets used for
functionality in the \code{signatureSearch} package, 
such as getting targets of query drugs in DrugBank, CLUE and STITCH databases, 
conducting GO enrichment analysis or computing WTCS p-values or tau scores 
of `gess_lincs` method.
}
\details{
The datasets in this package are:
  \itemize{
    \item \code{\link{cmap}} 
    \item \code{\link{cmap_expr}} 
    \item \code{\link{cmap_rank}} 
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
ah <- AnnotationHub()
## Load cmap database
# cmap <- ah[["AH69090"]]
## Load dtlink SQLite database
# conn <- load_sqlite("AH69083")
## Load RDS files
# goAnno <- ah[["AH69084"]]
}

\keyword{package}

