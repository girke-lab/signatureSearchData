\name{cmap_rank}
\alias{cmap_rank}
\title{Cmap Rank Signature Database}
\description{
CMap2 (Version build02) contains GESs for 1,309 drugs and eight cell lines that 
were generated with Affymetrix Gene Chips as expression platform. In some cases 
this includes drug treatments at different concentrations and time points. For 
consistency, the CMap2 data was reduced to drug treatments with concentrations 
and time points that are comparable to those used for the LINCS data. CMap2 data 
can be downloaded from GEO or its project site either in raw format or as rank 
transformed matrix. The ranks are based on DEG analyses of drug treatments 
(drug vs. no-drug) where the resulting Z-scores were used to generate the 
rank matrix. The latter was used here and is referred to as rankMatrix. The 
Affymetrix probe set identifiers stored in the row name slot of this matrix 
were translated into gene identifies. To obtain a matrix with unique gene 
identifiers, the ranks for genes represented by more than one probe set were 
averaged and then re-ranked accordingly. This final gene level rank matrix,
referred to as \code{cmap_rank}, contains rank profiles for 12,403 genes from 
1,309 compound treatments in up to 5 cells corresponding to a total of 
3,587 treatment signatures. This matrix can be used for all GESS methods in 
the signatureSearch package that are compatible with rank data, such as the 
\code{gess_cmap} method.
}
\details{
The \code{cmap_rank} data can be downloaded from Bioconductor’s AnnotationHub as 
HDF5 file. Since CMap2 is much smaller than LINCS, it can be imported in its 
entirety into a SummarizedExperiment object without excessive memory requirements
as shown in the Examples section.

The loaded \code{cmap_rank} data object is generated from the \code{rankMatrix}
downloaded from the CMap project site.
For documentation and code of generating the cmap_rank databases from sources, 
please refer to the vignette of this package by running 
\code{browseVignettes("signatureSearchData")} in R.
}

\references{
CMap project site: https://portals.broadinstitute.org/cmap
}

\examples{
library(AnnotationHub)
#ah <- AnnotationHub()
#query(ah, c("signatureSearchData", "cmap_rank"))
# cmap_rank_path <- ah[["AH69763"]]
# se <- readHDF5chunk(h5file=cmap_rank_path, colindex=1:3587)
}

\keyword{datasets}

