#' Calculate Mean Expression Values of LINCS Level 3 Data
#' 
#' Function calculates mean expression values for replicated samples of LINCS
#' Level 3 data that have been treated by the same compound in the same cell type
#' at a chosen concentration and treatment time. Usually, the function is used
#' after filtering the Level 3 data with \code{inst_filter}. The results (here
#' matrix with mean expression values) are saved to an HDF5 file. The latter is
#' referred to as the `lincs_expr` database.
#' @param gctx character(1), path to the LINCS Level 3 gctx file
#' @param inst tibble, LINCS Level 3 instances after filtering for specific concentrations and times
#' @param h5file character(1), path to the destination HDF5 file
#' @param chunksize number of columns of the matrix to be processed at a time to limit memory usage
#' @param overwrite TRUE or FALSE, whether to overwrite or append data to an
#' existing 'h5file'
#' @return HDF5 file, representing the 'lincs_expr' database
#' @importFrom signatureSearch createEmptyH5
#' @importFrom signatureSearch parse_gctx
#' @importFrom signatureSearch append2H5
#' @importFrom rhdf5 h5ls
#' @import AnnotationHub
#' @examples
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", package="signatureSearchData")
#' h5file <- tempfile(fileext=".h5")
#' inst <- data.frame(inst_id=c("ASG001_MCF7_24H:BRD-A79768653-001-01-3:10",
#'                              "CPC012_SKB_24H:BRD-K81418486:10"), 
#'     pert_cell_factor=c('sirolimus__MCF7__trt_cp', 'vorinostat__SKB__trt_cp'))
#' meanExpr2h5(gctx, inst, h5file, overwrite=TRUE)
#' @export
#' 
meanExpr2h5 <- function(gctx, inst, h5file, chunksize=2000, overwrite=TRUE){
    ## Get mean expression values of replicate samples from the same drug treatment in cells
    pert_cell_list <- split(inst$inst_id, inst$pert_cell_factor)
    chunk_list <- suppressWarnings(
        split(pert_cell_list, rep(seq_len(ceiling(length(pert_cell_list)/chunksize)), each=chunksize)))
    ## Creat an empty h5file
    if(file.exists(h5file)){
        if(overwrite){
            createEmptyH5(h5file, delete_existing=TRUE)
        }
    } else {
        createEmptyH5(h5file, delete_existing=FALSE)
    }
    ## append mean expr mat in each chunk to h5file
    lapply(chunk_list, function(pcl){
        cid_all <- unlist(pcl)
        mat_all <- parse_gctx(gctx, cid=cid_all, matrix_only=TRUE)
        mat_all <- mat_all@mat
        expr <- vapply(pcl, function(cid, mat_all){
            rowMeans(as.matrix(mat_all[,cid]))
        }, mat_all=mat_all, FUN.VALUE=numeric(12328))
        append2H5(x=expr, h5file, printstatus=FALSE)
    })
    h5ls(h5file)
}