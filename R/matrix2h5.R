#' Write Matrix to HDF5 file
#' 
#' Function writes matrix object to an HDF5 file.
#' @param matrix matrix to be written to HDF5 file, row and column name slots need
#' to be populated
#' @param h5file character(1), path to the hdf5 destination file
#' @param overwrite TRUE or FALSE, whether to overwrite or append
#' matrix to an existing 'h5file'
#' @return HDF5 file containing exported matrix
#' @examples 
#' mat <- matrix(rnorm(12), nrow=3, dimnames=list(paste0("r",1:3), paste0("c",1:4)))
#' h5file <- tempfile(fileext=".h5")
#' matrix2h5(matrix=mat, h5file=h5file, overwrite=TRUE)
#' @export
matrix2h5 <- function(matrix, h5file, overwrite=TRUE){
    if(file.exists(h5file)){
        if(isTRUE(overwrite)){
            createEmptyH5(h5file, delete_existing=TRUE)
        }
    } else {
        createEmptyH5(h5file, delete_existing=FALSE)
    }
    append2H5(matrix, h5file)
}
