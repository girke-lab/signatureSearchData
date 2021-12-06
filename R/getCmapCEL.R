#' This function will download the 7,056 CEL files from the CMap2 project 
#' site (http://www.broadinstitute.org/cmap), and save each of them to a
#' subdirectory named CEL under data. Since this download step will take 
#' a long time, the \code{rerun} argument has been assigned \code{FALSE} in
#' the example code below to avoid running it accidentally. If the raw data are not
#' needed, users can skip this time consuming step and work with the preprocessed
#' \code{\link{cmap}} or \code{\link{cmap_expr}} database downloaded from the 
#' ExperimentHub instead.
#' @title Download CMap2 CEL Files
#' @param rerun TRUE or FALSE, whether to download the data
#' @param dest_dir character(1), path to the destination directory
#' @return download files
#' @examples 
#' getCmapCEL(rerun=FALSE) # set 'rerun' to TRUE if download
#' @importFrom R.utils bunzip2
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @export
getCmapCEL <- function(rerun=TRUE, dest_dir=tempdir()) {
    if(rerun) {
        ## Check for presence of target directory
        if(!dir.exists(dest_dir))
            stop(paste0("Target directory ", dest_dir, " does not exist. 
             Create it to run function."))
        dest_cel <- file.path(dest_dir, "CEL")
        dir.create(dest_cel)
        ## Download CEL files
        for(vol in 1:7){
            dest_vol <- paste0(dest_cel, "/cmap_build02.volume", vol, "of7.zip")
            download.file(paste0("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume", vol, "of7.zip"), 
                          dest_vol)
            unzip(dest_vol, exdir=dest_cel)
            unlink(dest_vol)
        }
        
        ## Uncompress CEL files
        myfiles <- list.files(dest_cel, pattern=".CEL.bz2$", full.names=TRUE)
        for(i in myfiles) R.utils::bunzip2(i, remove=TRUE)
    } else {
        message("To execute function, set 'rerun=TRUE'")
    }
}