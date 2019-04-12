## Create empty HDF5 file
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

## Write in chunks to existing/empty HDF5 file
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
#' Read in a subset of matrix from HDF5 file as 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#' 
#' @param h5file character(1), path to the HDF5 file
#' @param colindex index of columns of the matrix to be read in
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", package="signatureSearchData")
#' h5file <- tempfile(fileext=".h5")
#' gctx2h5(gctx, cid=1:2, new_cid=c('sirolimus__MCF7__trt_cp', 'vorinostat__SKB__trt_cp'), 
#'         h5file=h5file, chunksize=5000, overwrite=TRUE)
#' se <- readHDF5chunk(h5file, colindex=1:2)
#' @export
#' 
readHDF5chunk <- function(h5file, colindex=seq_len(10)) {
    m <- h5read(h5file, "assay", index=list(NULL, colindex))
    mycol <- h5read(h5file, "colnames", index=list(colindex, 1))
    myrow <- h5read(h5file, "rownames")
    rownames(m) <- as.character(myrow[,1])
    colnames(m) <- as.character(mycol[,1])
    se <- SummarizedExperiment(assays=list(score=m))
    return(se)
}





