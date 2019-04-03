#' Load "dtlink_db_clue_sti" SQLite database into R by returning the
#' connection to the SQLite file.
#' 
#' The dtlink database is stored in \code{\link{AnnotationHub}}
#' @title Load dtlink database
#' @param ah_id AnnotationHub unique identifiers, of the form "AH12345",
#' of the hub records.
#' @return SQLiteConnection 
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom RSQLite SQLite
#' @importFrom RSQLite dbConnect
#' @examples
#' library(AnnotationHub)
#' #load_sqlite("AH69083")
#' @seealso 
#' \code{\link{dtlink_db_clue_sti}}
#' @export
load_sqlite <- function(ah_id){
    ah <- AnnotationHub()
    path <- ah[[ah_id]]
    #path <- fileName(ah[ah_id])
    conn <- dbConnect(SQLite(), path)
    return(conn)
}

#' Load signature databases (cmap etc.), which are HDF5 backed 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object and stored
#' in \code{\link{AnnotationHub}}
#' 
#' @title Load signature database
#' @param ah_id_h5 AnnotationHub id representing 'assays.h5' file from 
#' \code{\link[HDF5Array]{saveHDF5SummarizedExperiment}} function.
#' @param ah_id_rds AnnotationHub id representing 'se.rds' file from 
#' \code{\link[HDF5Array]{saveHDF5SummarizedExperiment}} function.
#' @param db_name character(1), name of the signatur database, e.g., "cmap"
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#' with `assays` slot as \code{\link[DelayedArray]{DelayedMatrix}}
#' @seealso 
#' \code{\link{cmap}}, \code{\link{lincs}}
#' @importFrom AnnotationHub hubCache
#' @import HDF5Array
#' @examples 
#' library(AnnotationHub)
#' #se <- load_sigdb("AH69075","AH69076","cmap")
#' @export
load_sigdb <- function(ah_id_h5, ah_id_rds, db_name){
    ah <- AnnotationHub()
    h5_path = ah[[ah_id_h5]]
    rds_path = ah[[ah_id_rds]]
    hub_dir <- hubCache(ah)
    wd <- getwd()
    setwd(hub_dir)    
    # Create db directory and symbolic link to the h5, rds file under this dir
    if(!dir.exists(db_name)) dir.create(db_name)
    suppressWarnings(file.symlink(c(h5_path, rds_path), 
                 file.path(db_name, c("assays.h5", "se.rds"))))
    se <- loadHDF5SummarizedExperiment(db_name)
    setwd(wd)
    return(se)
}

