#' Create empty HDF5 file 
#' 
#' This function can be used to create an empty HDF5 file by defining the file
#' path and compression level. The created HDF5 file has three datasets named 
#' 'assay', 'colnames', 'rownames' under root group. 
#' 'assay' is used to store values a numeric matrix, 'colnames' and 'rownames' 
#' are used to store character vectors of column names and row names of the 
#' matrix, respectively.
#' 
#' @param h5file character(1), path to the HDF5 file to be created
#' @param delete_existing logical, whether to delete the existing HDF5 file
#' @param level The compression level used. An integer value between 0 
#' (no compression) and 9 (highest and slowest compression).
#' @return empty HDF5 file
#' @examples
#' tmp_file <- tempfile(fileext=".h5")
#' createEmptyH5(tmp_file, level=6)
#' @export
createEmptyH5 <- function(h5file, delete_existing=FALSE, level=6) {
    if(delete_existing==TRUE) unlink(h5file)
    h5createFile(file=h5file)
    h5createDataset(h5file, "assay", c(0,0), c(H5Sunlimited(), H5Sunlimited()), 
                    chunk=c(12328,1), level=level)
    h5createDataset(h5file, "colnames", c(0,1), c(H5Sunlimited(), 1), 
                    storage.mode='character', size=1000, level=level)
    h5createDataset(h5file, "rownames", c(0,1), c(H5Sunlimited(), 1), 
                    storage.mode='character', size=100, level=level)
}

#' Append matrix to an existing/empty HDF5 file
#' @param x matrix object in R to be appended. If the HDF5 file is not empty,
#' it needs to have the same row number as the matrix already existing in the HDF5
#' file, the columns will be appended to the existing columns. 
#' @param h5file character(1), path to the existing/empty HDF5 file
#' @param printstatus logical, whether to print status
#' @return existing/empty HDF5 file appended by a new matrix
#' @examples 
#' mat <- matrix(1:12, nrow=3)
#' rownames(mat) <- paste0("r", 1:3); colnames(mat) <- paste0("c", 1:4)
#' tmp_file <- tempfile(fileext=".h5")
#' createEmptyH5(tmp_file)
#' append2H5(mat, tmp_file)
#' rhdf5::h5ls(tmp_file)
#' @export
append2H5 <- function(x, h5file, printstatus=TRUE) {
    status <- h5ls(h5file)[c("name", "dim")]
    rowstatus <- as.numeric(gsub(" x \\d{1,}$", "", status[status$name=="assay", "dim"]))
    colstatus <- as.numeric(gsub("^\\d{1,} x ", "", status[status$name=="assay", "dim"]))
    nrows <- nrow(x) 
    ncols <- colstatus + ncol(x)
    h5set_extent(h5file, "assay", c(nrows, ncols))
    h5write(x, h5file, "assay", index=list(seq_len(nrows), (colstatus+1):ncols))
    h5set_extent(h5file, "colnames", c(ncols,1))
    h5write(colnames(x), h5file, "colnames", index=list((colstatus+1):ncols, 1))
    if(any(duplicated(h5read(h5file, "colnames")[,1]))) 
        warning("Column names contain duplicates!")
    h5set_extent(h5file, "rownames", c(nrows,1))
    h5write(rownames(x), h5file, "rownames", index=list(seq_len(nrows), 1))
    if(any(duplicated(h5read(h5file, "rownames")[,1]))) 
        warning("Row names contain duplicates!")
    if(printstatus==TRUE) h5ls(h5file, all=TRUE)[c("dim", "maxdim")]
    h5closeAll()
}

#' Read HDF5 file as SummarizedExperiment object
#' 
#' Read in a subset of matrix by setting column index from an HDF5 file as a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. The
#' HDF5 file need to has three datasets named 
#' 'assay', 'colnames', 'rownames' under root group. 
#' 'assay' is used to store values a numeric matrix, 'colnames' and 'rownames' 
#' are used to store character vectors of column names and row names of the 
#' matrix, respectively.
#' 
#' @param h5file character(1), path to the HDF5 file
#' @param colindex integer vector, index of the columns of the matrix to be read in
#' @param colnames character vector, names of the columns of the matrix to be 
#' read in. If 'colnames' is set, 'colindex' will be ignored.
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @seealso 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @examples
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", package="signatureSearchData")
#' h5file <- tempfile(fileext=".h5")
#' gctx2h5(gctx, cid=1:2, new_cid=c('sirolimus__MCF7__trt_cp', 'vorinostat__SKB__trt_cp'), 
#'         h5file=h5file, chunksize=5000, overwrite=TRUE)
#' se <- readHDF5chunk(h5file, colindex=1:2)
#' @export
#' 
readHDF5chunk <- function(h5file, colindex=seq_len(10), colnames=NULL) {
    if(! is.null(colnames)){
        all_trts <- h5read(h5file, "colnames", drop=TRUE)
        colindex2 <- which(all_trts %in% colnames)
        m <- h5read(h5file, "assay", index=list(NULL, colindex2))
        colindex <- colindex2
    } else {
        m <- h5read(h5file, "assay", index=list(NULL, colindex))
    }
    mycol <- h5read(h5file, "colnames", index=list(colindex, 1))
    myrow <- h5read(h5file, "rownames")
    rownames(m) <- as.character(myrow[,1])
    colnames(m) <- as.character(mycol[,1])
    if(! is.null(colnames)){
        m = m[,colnames, drop=FALSE]
    }
    se <- SummarizedExperiment(assays=list(score=m))
    h5closeAll()
    return(se)
}





