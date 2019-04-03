#' Read in the matrix stored in gctx file by chunks and save it as hdf5 backed 
#' 'DelayedArray'. 
#' 
#' The large gctx file was read in by chunks. It avoids consuming all of the
#' memories by reading the full data at once. The chunk matrix was stored as
#' \code{\link[HDF5Array]{HDF5Array}}, which was realized by hdf5 backend. 
#' Then the chunks were
#' combined as a \code{\link[DelayedArray]{DelayedMatrix}} object that contains
#'  all the data read from gctx file 
#' @title read large gctx file by chunks
#' @param gctx path to the gctx file downloaded and decompressed from GEO  
#' @param cid all the column ids to be read in
#' @param chunk number of columns read in at once
#' @return DelayedMatrix object containing all the columns in 'cid' read from 'gctx'
#' @seealso 
#' \code{\link[DelayedArray]{DelayedMatrix}}, \code{\link[HDF5Array]{HDF5Array}}
#' @examples 
#' a = 1
#' # lincs42_dlmat <- read_gctx(gctx="./data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", 
#' # cid=meta42$sig_id))
#' @export 
read_gctx <- function(gctx, cid, chunk=5000){
    cid_list <- suppressWarnings(
        split(cid, rep(seq_len(ceiling(length(cid)/chunk)), each=chunk)))
    delay_list <- lapply(cid_list, function(x){
        mat1 <- parse_gctx(gctx, cid=x, matrix_only=TRUE)
        mat1 <- mat1@mat
        hdf <- as(mat1, "HDF5Array")
        rownames(hdf) <- rownames(mat1)
        return(hdf)
    })
    delay <- do.call("cbind", delay_list)
    return(delay)
}