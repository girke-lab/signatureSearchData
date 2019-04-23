\name{cmap_rank}
\alias{cmap_rank}
\title{cmap rank signature database}
\description{
  The `cmap_rank` database represents the ranks of 12,403 genes after 
  differential expression analysis from 1,281 compound treatments in 4 cells 
  (3,478 signatures in total). It is stored in an HDF5 file.
}
\details{
The `cmap_rank` signature database can be used to calculate
the connectivity scores to the query gene sets by `gess_cmap` method in the
`signatureSearch` package.

The Connectivity Map \url{https://portals.broadinstitute.org/cmap/} group has 
pre-processed a `rankMatrix` storing the ranks of probe sets 
for each instance (a treatment and control pair and the probe sets were ordered 
by their extent of differential expression between this treatment and control pair).
The `rankMatrix` was transferred from probe set level to gene level, selected,
and stored in an HDF5 file (`cmap_rank` database).
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("signatureSearchData", "cmap_rank"))
# cmap_rank <- ah[["AH69763"]]
}

\keyword{datasets}

