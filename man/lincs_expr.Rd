\name{lincs_expr}
\alias{lincs_expr}
\title{LINCS Intensity Signature Database} 
\description{
The LINCS Level 3 data can be downloaded from GEO the same way as described
for the Level 5 data available at help file of \code{\link{lincs}}. 
The Level 3 data contain normalized gene expression values across all treatments 
and cell lines used by LINCS. The Level 3 signatures were filtered using the 
same dosage and duration criteria as the Level 5 data. The biological replicate 
information included in the Level 3 data were collapsed to mean values. 
Subsequently, the resulting matrix of mean expression values was written to an 
HDF5 file. The latter is referred to as \code{lincs_expr} database containing 
38,824 signatures for a total of 5,925 small molecule treatments and 30 cell lines. 
Although the LINCS Level 3 and 5 data are filtered here the same way, the 
number of small molecules represented in the Level 3 data (5,925) is smaller 
than in the Level 5 data (8,140). The reason for this inconsistency is most 
likely that the Level 3 dataset, downloadable from GEO, is incomplete.
}
\details{
The filtered and processed LINCS Level3 data (lincs_expr) can be loaded from 
Bioconductor’s AnnotationHub interface as shown in Examples section. 
In this case the loaded lincs_expr instance includes mean expression values of
12,328 genes for 5,925 compound treatments across a total of 30 cell lines
(38,824 signatures in total). 
This data set can be used by all correlation-based GESS methods provided by the
signatureSearch package.

The loaded \code{lincs_expr} data object is filtered and generated from the original
LINCS level 3 data stored at GEO
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742}. 
For documentation and code of generating the lincs_expr databases from sources, 
please refer to the vignette of this package by running 
\code{browseVignettes("signatureSearchData")} in R.
}

\references{
  LINCS source data at 5 levels from GEO: 
  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742} 
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("signatureSearchData", "lincs_expr"))
# lincs_expr_path <- ah[["AH69093"]]
}

\keyword{datasets}

