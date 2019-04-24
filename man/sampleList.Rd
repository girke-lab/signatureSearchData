% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleList.R
\name{sampleList}
\alias{sampleList}
\title{Generate drug treatment vs. control list}
\usage{
sampleList(cmap, myby)
}
\arguments{
\item{cmap}{cmap instances annotation}

\item{myby}{"CMP" or "CMP_CELL", "CMP": by compound treatments in all cells;
"CMP_CELL": by compound treatments in individual cell}
}
\value{
list
}
\description{
Generate CEL file list for treatment vs. control comparisons
}
\examples{
path <- system.file("extdata", "cmap_instances_02.txt", package="signatureSearchData")
cmap_inst <- read.delim(path, check.names=FALSE) 
comp_list <- sampleList(cmap_inst, myby="CMP_CELL")
}