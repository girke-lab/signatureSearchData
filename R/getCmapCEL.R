#' Download CEL files from CMap project \url{http://www.broadinstitute.org/cmap}
#' @title get CEL files
#' @param rerun TRUE or FALSE, whether to download the data
#' @return files
#' @examples 
#' getCmapCEL(rerun=FALSE) # set 'rerun' to TRUE if download
#' @importFrom R.utils bunzip2
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @export
getCmapCEL <- function(rerun=TRUE) {
    if(isTRUE(rerun)) {
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