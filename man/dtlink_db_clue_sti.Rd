\name{dtlink_db_clue_sti}
\alias{dtlink_db_clue_sti}
\title{Annotation of drug-target links from DrugBank, CLUE and STITCH databases} 
\description{
It is a SQLite database containing drug-target links in DrugBank, CLUE, 
and STITCH databases. 
}
\details{
It is an intermediate database used for functionalities in the 
\code{signatureSearch}
software package to get targets of query drugs in DrugBank 
\url{https://www.drugbank.ca/}, CLUE \url{https://clue.io/} and STITCH
\url{http://stitch.embl.de/} databases.

The STITCH database provides confidence scores for each drug-target interaction. 
Only interactions with confidence scores >=0.7 were selected. 
To minimize noise, e.g. from promiscuous binders, drugs with more than 100 
distinct targets were removed. Similarly, targets with more than 100 annotated 
drugs were excluded.
}

\seealso{
  \itemize{
    \item \code{\link{signatureSearch}}
    \item \code{\link[signatureSearch]{get_targets}}
  }
}

\keyword{datasets}

