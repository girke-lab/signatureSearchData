\name{lincs}
\alias{lincs}
\title{LINCS Z-score Signature Database} 
\description{
The L1000 assay, used for generating the LINCS data, measures the expression of 
978 landmark genes and 80 control genes by loading amplified mRNA populations 
onto beads and then detecting their abundance with a fluorescent-based method. 
The expression of 11,350 additional genes is imputed from the landmark genes 
by using as training data a large collection of Affymetrix gene chips.

The LINCS data have been pre-processed by the Broad Institute to 5 different 
levels and are available for download from GEO. Level 1 data are the raw mean 
fluorescent intensity values that come directly from the Luminex scanner. 
Level 2 data are the expression intensities of the 978 landmark genes. 
They have been normalized and used to impute the expression of the additional 
11,350 genes, forming Level 3 data. A robust z-scoring procedure was used to 
generate differential expression values from the normalized profiles (Level 4). 
Finally, a moderated z-scoring procedure was applied to the replicated samples 
of each experiment (mostly 3 replicates) to compute a weighted average signature 
(Level 5). For a more detailed description of the preprocessing methods used by 
the LINCS project, readers want to refer to the LINCS user guide.

Disregarding replicates, the LINCS data set contains 473,647 signatures with 
unique cell type and treatment combinations. This includes 19,811 drug-like 
small molecules tested on different cell lines at multiple concentrations and 
treatment times. In addition to compounds, several thousand genetic perturbations 
(gene knock-downs and over expressions) have been tested. Currently, the data 
stored in this package is restricted to signatures of small molecule treatments 
across different cells lines. However, users have the option to assemble any 
custom collection of the LINCS data. For consistency, only signatures at one 
specific concentration (10uM) and one time point (24h) have been selected for 
each small molecule in the default collection. These choices are similar to the 
conditions used in primary high-throughput compound screens of cell lines. 
Since the selected compound concentrations and treatment duration have not been 
tested by LINCS across all cell types yet, a subset of compounds had to be 
selected that best met the chosen treatment requirements. This left us with 
8,104 compounds that were uniformly tested at the chosen concentration and 
treatment time, but across variable numbers of cell lines. The total number of
expression signatures meeting this requirement is 45,956, while the total number 
of cell lines included in this data set is 30.
}
\details{
The LINCS sub-dataset, filtered and assembled according to the above criteria, 
can be downloaded from Bioconductor’s \code{ExperimentHub} as HDF5 file. 
In this case the loaded data instance includes moderated Z-scores from 
DE analyses of 12,328 genes for 8,140 compound treatments across a total of 
30 cell lines corresponding to 45,956 expression signatures. This data set can 
be used by all set-based and correlation-based GESS methods provided by the 
signatureSearch package.

The loaded \code{lincs} data object is filtered and generated from the original
LINCS level 5 data stored at GEO:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742. 
For documentation and code of generating the lincs databases from sources, 
please refer to the vignette of this package by running 
\code{browseVignettes("signatureSearchData")} in R.
}

\references{
  LINCS source data at 5 levels from GEO: 
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
  
  LINCS User Guide v2.1: 
  https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#
}

\examples{
library(ExperimentHub)
# eh <- ExperimentHub()
# query(eh, c("signatureSearchData", "lincs"))
# lincs_path <- eh[['EH3226']]
# rhdf5::h5ls(lincs_path) 
}

\keyword{datasets}

