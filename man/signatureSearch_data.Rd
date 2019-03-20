\name{signatureSearch_data}
\alias{signatureSearch_data}
\title{Annotation data package for signatureSearch}
\description{
This is an annotation data package for \code{signatureSearch} softeware packge. 
It stores the CMAP/LINCS reference databases used for GESS methods in the
\code{signatureSearch} package as well as other intermediate datasets used for
functionalities in the \code{signatureSearch} package, 
such as getting targets of query drugs in DrugBank, CLUE and STITCH databases, 
conducting GO enrichment analysis or computing p-values or Tau scores 
in the results from `gess_lincs` method.
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
For codes and details of generating the datasets, please refer to the 
\code{make-data.R} script at 
\url{https://github.com/yduan004/signatureSearch_data/blob/master/inst/scripts/make-data.R}
}

\seealso{
   \code{\link{signatureSearch}}
}

\keyword{package}

