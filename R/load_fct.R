#' This function loads the \code{\link{dtlink_db_clue_sti}} SQLite database 
#' into R by downloading the database from AnnotationHub if it is not in the cache
#' file and creating the connection to the SQLite file in R via \code{dbConnect} 
#' function in \code{RSQLite} package.
#' 
#' @title Load SQLite Database from AnnotationHub
#' @param ah_id AnnotationHub unique identifiers, of the form "AH12345",
#' of the hub records.
#' @return SQLiteConnection 
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom RSQLite SQLite
#' @importFrom RSQLite dbConnect
#' @examples
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' qr <-query(ah, c("signatureSearchData", "dtlink"))
#' conn <- load_sqlite("AH69083")
#' RSQLite::dbListTables(conn)
#' RSQLite::dbDisconnect(conn)
#' @seealso 
#' \code{\link{dtlink_db_clue_sti}}, \code{\link[AnnotationHub]{AnnotationHub}}
#' @export
load_sqlite <- function(ah_id){
    ah <- AnnotationHub()
    path <- suppressMessages(ah[[ah_id]])
    conn <- dbConnect(SQLite(), path)
    return(conn)
}

