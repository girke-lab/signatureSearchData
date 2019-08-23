\name{cmap}
\alias{cmap}
\title{CMap2 LFC Signature Database}
\description{
The MAS5 normalized CEL files from the help file of \code{\link{cmap_expr}} 
can be used for differential expression (DE) analysis with `limma` package to get a matrix 
containing the LFC (log2 fold change) scores. The columns of the matrix
are the treatment v.s. control instances, which is defined in the 
'cmap_instances_02.txt' file downloaded from CMap project website. 
The same as the `cmap_rank` database,
only one treatment condition is selected for a compound in a cell. 
The resulting matrix represents LFC of 12,403 genes for 1,281 compound treatments 
in up to 5 cells (3,478 signatures in total). The latter was stored in an HDF5 
file, which is referred to as the \code{cmap} database. Note, The number of 
compound treatments in \code{cmap} database is slightly different from that of 
the \code{cmap_expr} database. The reason is that some of the compound treatment 
is discarded if the number of control and treatment samples are less than 3 during the DE analysis.
}
\details{
The cmap database can be downloaded as HDF5 file from Bioconductor’s 
AnnotationHub as shown in the Example section. 

The same as the \code{cmap_expr} dataset, the loaded \code{cmap} data object is 
also generated from the raw CEL files downloaded from the CMap project site and
processed as describe above.
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("signatureSearchData", "cmap"))
# cmap_path <- ah[["AH69090"]]
}

\keyword{datasets}

