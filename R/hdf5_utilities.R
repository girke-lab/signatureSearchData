#' Create Empty HDF5 File 
#' 
#' This function can be used to create an empty HDF5 file where the user defines
#' the file path and compression level. The empty HDF5 file has under its root
#' group three data slots named 'assay', 'colnames' and 'rownames' for storing a
#' \code{numeric matrix} along with its column names (\code{character}) and row names
#' (\code{character}), respectively.
#' 
#' @param h5file character(1), path to the HDF5 file to be created
#' @param delete_existing logical, whether to delete an existing HDF5 file with identical path
#' @param level The compression level used, here given as integer value between 0
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

#' Append Matrix to HDF5 File
#' 
#' Function to write matrix data to an existing HDF5 file. If the file contains 
#' already matrix data then both need to have the same number of rows. The append
#' will be column-wise.
#' @param x matrix object to write to an HDF5 file. If the HDF5 file is not empty,
#' the exported matrix data needs to have the same number rows as the matrix
#' stored in the HDF5 file, and will be appended column-wise to the existing
#' one.
#' @param h5file character(1), path to existing HDF5 file that can be empty or
#' contain matrix data
#' @param printstatus logical, whether to print status
#' @return HDF5 file storing exported matrix
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

#' Import HDF5 Data into SummarizedExperiment Object
#' 
#' Imports user-definable subsets of matrix data from an HDF5 file into a 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. The
#' corresponding HDF5 file is expected to have three data components named
#' 'assay', 'colnames' and 'rownames' containing the numeric values, column
#' names and row names of a matrix, respectively.
#' 
#' @param h5file character(1), path to HDF5 file
#' @param colindex integer vector, position index of the matrix columns to be imported
#' @param colnames character vector, names of the columns of the matrix to be 
#' imported. If 'colnames' is set, 'colindex' will be ignored.
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





