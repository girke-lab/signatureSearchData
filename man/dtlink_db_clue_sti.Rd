\name{dtlink_db_clue_sti}
\alias{dtlink_db_clue_sti}
\title{Drug-Target Annotation Resources}
\description{
This is a SQLite database containing drug-target links (drug names to gene SYMBOL ids)
obtained from DrugBank, CLUE, and STITCH annotation databases. 
}
\details{
It is an intermediate database used for functionality in the 
\code{signatureSearch} software package to get target genes/proteins of query 
drugs in DrugBank (https://www.drugbank.ca/), CLUE (https://clue.io/) and STITCH
(http://stitch.embl.de/) databases.

The target proteins of drugs are mapped to their encoding genes by using 
identifier mapping resources from R/Bioconductor such as the \code{org.Hs.eg.db}
annotation package. Thus, the drug-target links in this database are drug names
to gene SYMBOLs mapping. 

The STITCH database provides confidence scores for each drug-target interaction.
Only interactions with confidence scores >=0.7 were selected.

The drug-target interaction annotations from the above three databases are combined.
To minimize noise, e.g. from promiscuous binders, drugs with more than 100 
distinct targets are removed. Similarly, targets with more than 100 annotated 
drugs are excluded.
}

\seealso{
    \code{get_targets}
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "dtlink_db_clue_sti"))
names(qr)
}

\keyword{datasets}

