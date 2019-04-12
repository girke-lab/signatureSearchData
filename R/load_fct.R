#' Load "dtlink_db_clue_sti" SQLite database into R by returning the
#' connection to the SQLite file.
#' 
#' The `dtlink_db_clue_sti` database is stored in \code{\link{AnnotationHub}}
#' @title Load dtlink database
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

