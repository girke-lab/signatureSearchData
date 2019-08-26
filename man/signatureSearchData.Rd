\name{signatureSearchData}
\alias{signatureSearchData}
\title{signatureSearchData: Reference Data for Gene Expression Signature Searching}
\description{
This package provides access to the reference data used by the associated 
signatureSearch software package. The latter allows to search with a query 
gene expression signature (GES) a database of treatment GESs to identify 
cellular states sharing similar expression responses (connections). This way 
one can identify drugs or gene knockouts that induce expression phenotypes 
similar to a sample of interest. The resulting associations may lead to novel 
functional insights how perturbagens of interest interact with biological systems.

Currently, this data package includes GES data from the CMap (Connectivity Map) 
and LINCS (Library of Network-Based Cellular Signatures) projects that are 
largely based on drug and genetic perturbation experiments performed on 
variable numbers of human cell lines (Lamb et al. 2006; Subramanian et al. 2017).
In signatureSearchData these data sets have been preprocessed to be compatible 
with the different gene expression signature search (GESS) algorithms 
implemented in signatureSearch. The preprocessed data types include but are not 
limited to normalized gene expression values (e.g. intensity values), log fold 
changes (LFC) and Z-scores, p-values or FDRs of differentially expressed genes 
(DEGs), rankings based on selected preprocessing routines or sets of top 
up/down-regulated DEGs.

The CMap data were downloaded from the CMap project site (Version build02). 
The latter is a collection of over 7,000 gene expression profiles (signatures) 
obtained from perturbation experiments with 1,309 drug-like small molecules 
on five human cancer cell lines. The Affymetrix Gene Chip technology was used 
to generate the CMAP2 data set.

In 2017, the LINCS Consortium generated a similar but much larger data set 
where the total number of gene expression signatures was scaled up to over 
one million. This was achieved by switching to a much more cost effective gene 
expression profiling technology called L1000 assay (Peck et al. 2006; Edgar, 
Domrachev, and Lash 2002). The current set of perturbations covered by the 
LINCS data set includes 19,811 drug-like small molecules applied at variable 
concentrations and treatment times to ~70 human non-cancer (normal) and cancer 
cell lines. Additionally, it includes several thousand genetic perturbagens 
composed of gene knockdown and over-expression experiments.

The data structures and search algorithms used by signatureSearch and 
signatureSearchData are designed to work with most genome-wide expression data 
including hybridization-based methods, such as Affymetrix or L1000, as well as 
sequencing-based methods, such as RNA-Seq. Currently, signatureSearchData does 
not include preconfigured RNA-Seq reference data mainly due to the lack of 
large-scale perturbation studies (e.g. drug-based) available in the public 
domain that are based on RNA-Seq. This situation may change in the near future 
once the technology has become more affordable for this purpose.

This package also contains other intermediate annotation datasets, such as drug-target 
information of small molecules in DrugBank, CLUE and STITCH databases used for 
getting target annotation of query drugs, GO system annotations and drugs to GO 
functional category mappings used for signatureSearch's functional enrichment 
analysis (FEA) routines or null distribution of Enrichment 
Scores (ES) used for computing WTCS p-values from LINCS GESS method.
}
\details{
This package contains the following main datasets:
  \itemize{
    \item \code{\link{cmap}} 
    \item \code{\link{cmap_expr}} 
    \item \code{\link{cmap_rank}} 
    \item \code{\link{lincs}} 
    \item \code{\link{lincs_expr}} 
    \item \code{\link{dtlink_db_clue_sti}}
    \item \code{\link{taurefList}}
  }
For documentation of generating the CMAP/LINCS databases from sources, 
please refer to the vignette of this package by running 
\code{browseVignettes("signatureSearchData")} in R.

This package also includes several functions used for creating an empty HDF5 file,
appending a matrix object in batches to an HDF5 file, as well as read in an HDF5 file
into a \code{SummarizedExperiment} object by selecting columns.

The main functions include:
\itemize{
    \item \code{\link{createEmptyH5}} 
    \item \code{\link{append2H5}} 
    \item \code{\link{readHDF5chunk}} 
  }
}

\seealso{
   \code{signatureSearch}
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
ssd <- query(ah, c("signatureSearchData"))
ssd
ssd$title
as.list(ssd)[10]

## Load CMAP databases
# cmap_path <- ah[["AH69090"]] # stores LFC scores
# cmap_expr_path <- ah[["AH69091"]] # stores gene expression intensity values
# cmap_rank_path <- ah[["AH69763"]] # stores gene ranks
# se <- readHDF5chunk(h5file=cmap_rank_path, colindex=1:3587)
# 
## Load LINCS databases
# lincs_path <- ah[['AH69092']] # stores moderated z-scores
# h5ls(lincs_path)
# lincs_expr_path <- ah[['AH69093']] # stores gene expression intensity values
# 
## Load dtlink SQLite database
# conn <- load_sqlite("AH69083")
# RSQLite::dbListTables(conn)
# 
## Load GO annotation system (RDS files)
# goAnno <- ah[["AH69084"]]
}

\references{
Lamb, Justin, Emily D Crawford, David Peck, Joshua W Modell, Irene C Blat, 
Matthew J Wrobel, Jim Lerner, et al. 2006. “The Connectivity Map: Using 
Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease.” 
Science 313 (5795): 1929–35. http://dx.doi.org/10.1126/science.1132939.

Subramanian, Aravind, Rajiv Narayan, Steven M Corsello, David D Peck, Ted E 
Natoli, Xiaodong Lu, Joshua Gould, et al. 2017. “A Next Generation Connectivity 
Map: L1000 Platform and the First 1,000,000 Profiles.” Cell 171 (6): 
1437–1452.e17. http://dx.doi.org/10.1016/j.cell.2017.10.049.

Edgar, Ron, Michael Domrachev, and Alex E Lash. 2002. “Gene Expression Omnibus: 
NCBI Gene Expression and Hybridization Array Data Repository.” Nucleic Acids 
Res. 30 (1): 207–10. https://www.ncbi.nlm.nih.gov/pubmed/11752295.

Peck, David, Emily D Crawford, Kenneth N Ross, Kimberly Stegmaier, Todd R Golub,
and Justin Lamb. 2006. “A Method for High-Throughput Gene Expression Signature 
Analysis.” Genome Biol. 7 (7): R61. http://dx.doi.org/10.1186/gb-2006-7-7-r61.
}
\keyword{package}

