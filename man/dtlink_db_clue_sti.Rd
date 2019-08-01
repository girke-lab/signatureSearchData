\name{dtlink_db_clue_sti}
\alias{dtlink_db_clue_sti}
\title{Annotation of drug-target links from DrugBank, CLUE and STITCH databases}
\description{
It is a SQLite database containing drug-target links in DrugBank, CLUE, 
and STITCH databases. 
}
\details{
It is an intermediate database used for functionality in the 
\code{signatureSearch} software package to get targets of query drugs in 
DrugBank \url{https://www.drugbank.ca/}, CLUE \url{https://clue.io/} and STITCH
\url{http://stitch.embl.de/} databases.

The STITCH database provides confidence scores for each drug-target interaction.
Only interactions with confidence scores >=0.7 were selected. 
To minimize noise, e.g. from promiscuous binders, drugs with more than 100 
distinct targets were removed. Similarly, targets with more than 100 annotated 
drugs were excluded.
}

\seealso{
    get_targets
}

\examples{
library(AnnotationHub)
ah <- AnnotationHub()
qr <- query(ah, c("signatureSearchData", "dtlink_db_clue_sti"))
names(qr)
conn <- load_sqlite("AH69083")
RSQLite::dbListTables(conn)
RSQLite::dbDisconnect(conn)
}

\keyword{datasets}

