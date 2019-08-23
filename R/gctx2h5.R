#' Read matrix-like data from large gctx file in chunks and write result back to an 
#' HDF5 file.
#' @title Convert GCTX to HDF5 File
#' @param gctx character(1), path to gctx file from LINCS
#' @param cid character or integer vector referencing the
#' columns of the matrix to include
#' @param new_cid character vector of the same length as cid, assigning new
#' column names to matrix
#' @param h5file character(1), path of the hdf5 destination file
#' @param chunksize number of columns to import in each iteration to limit memory usage
#' @param overwrite TRUE or FALSE, whether to overwrite or to append to existing 'h5file'
#' @return HDF5 file
#' @import rhdf5
#' @examples 
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", package="signatureSearchData")
#' h5file <- tempfile(fileext=".h5")
#' gctx2h5(gctx, cid=1:2, new_cid=c('sirolimus__MCF7__trt_cp', 'vorinostat__SKB__trt_cp'), 
#'         h5file=h5file, chunksize=5000, overwrite=TRUE)
#' @export
#' 
gctx2h5 <- function(gctx, cid, new_cid=cid, h5file, chunksize=5000, 
                    overwrite=TRUE){
    cid_list <- suppressWarnings(
        split(cid, rep(seq_len(ceiling(length(cid)/chunksize)), each=chunksize)))
    new_cid_list <- suppressWarnings(
        split(new_cid, rep(seq_len(ceiling(length(new_cid)/chunksize)), each=chunksize)))
    if(file.exists(h5file)){
        if(isTRUE(overwrite)){
            createEmptyH5(h5file, delete_existing=TRUE)
        }
    } else {
        createEmptyH5(h5file, delete_existing=FALSE)
    }
    lapply(seq_along(cid_list), function(i){
        mat <- parse_gctx(gctx, cid=cid_list[[i]], matrix_only=TRUE)
        mat <- mat@mat
        colnames(mat) <- new_cid_list[[i]]
        append2H5(x=mat, h5file, printstatus=FALSE)
    })    
    h5ls(h5file)
}