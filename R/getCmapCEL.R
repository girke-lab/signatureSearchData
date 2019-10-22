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
#' @return download files
#' @examples 
#' getCmapCEL(rerun=FALSE) # set 'rerun' to TRUE if download
#' @importFrom R.utils bunzip2
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @export
getCmapCEL <- function(rerun=TRUE) {
    if(rerun) {
        ## Check for presence of target directory
        if(!dir.exists("./data/CEL"))
            stop("Target directory ./data/CEL does not exist. 
             Create it to run function.") 
        ## Download CEL files
        download.file("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume1of7.zip", 
                      "./data/CEL/cmap_build02.volume1of7.zip")
        unzip("./data/CEL/cmap_build02.volume1of7.zip", exdir="./data/CEL")
        unlink("./data/CEL/cmap_build02.volume1of7.zip")
        download.file("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume2of7.zip", 
                      "./data/CEL/cmap_build02.volume2of7.zip")
        unzip("./data/CEL/cmap_build02.volume2of7.zip", exdir="./data/CEL")
        unlink("./data/CEL/cmap_build02.volume2of7.zip")
        download.file("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume3of7.zip", 
                      "./data/CEL/cmap_build02.volume3of7.zip")
        unzip("./data/CEL/cmap_build02.volume3of7.zip", exdir="./data/CEL")
        unlink("./data/CEL/cmap_build02.volume3of7.zip")
        download.file("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume4of7.zip", 
                      "./data/CEL/cmap_build02.volume4of7.zip")
        unzip("./data/CEL/cmap_build02.volume4of7.zip", exdir="./data/CEL")
        unlink("./data/CEL/cmap_build02.volume4of7.zip")
        download.file("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume5of7.zip", 
                      "./data/CEL/cmap_build02.volume5of7.zip")
        unzip("./data/CEL/cmap_build02.volume5of7.zip", exdir="./data/CEL")
        unlink("./data/CEL/cmap_build02.volume5of7.zip")
        download.file("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume6of7.zip", 
                      "./data/CEL/cmap_build02.volume6of7.zip")
        unzip("./data/CEL/cmap_build02.volume6of7.zip", exdir="./data/CEL")
        unlink("./data/CEL/cmap_build02.volume6of7.zip")
        download.file("ftp://ftp.broad.mit.edu/pub/cmap/cmap_build02.volume7of7.zip", 
                      "./data/CEL/cmap_build02.volume7of7.zip")
        unzip("./data/CEL/cmap_build02.volume7of7.zip", exdir="./data/CEL")
        unlink("./data/CEL/cmap_build02.volume7of7.zip")
        
        ## Uncompress CEL files
        myfiles <- list.files("./data/CEL", pattern=".CEL.bz2$", full.names=TRUE)
        for(i in myfiles) R.utils::bunzip2(i, remove=TRUE)
    } else {
        message("To execute function, set 'rerun=TRUE'")
    }
}