#' Calculate mean expression value of compound treatment samples and save to HDF5 file
#' 
#' Calculate mean expression values of replicated samples that are treated by the same
#' compound in the same cell type with selected concentation and treatment time after 
#' filtering the LINCS level3 instances/samples. 
#' Then write the mean expression matrix to an HDF5 file, which is the `lincs_expr`
#' database.
#' @param gctx character(1), path to the LINCS level3 gctx file
#' @param inst tibble, LINCS level3 instances after filtering at specific concentration and time
#' @param h5file character(1), path to the destination hdf5 file
#' @param chunksize size of the matrix (number of columns) to be saved in HDF5 file by chunks
#' @param overwrite TRUE or FALSE, whether to overwrite or append the
#' matrix to the exsisting 'h5file'
#' @return HDF5 file, representing the 'lincs_expr' database
#' @examples
#' \dontrun{
#' meanExpr2h5(gctx="./data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",
#' inst=inst42_filter, h5file="./data/lincs_expr.h5")
#' }
#' @export
#' 
meanExpr2h5 <- function(gctx, inst, h5file, chunksize=2000, overwrite=TRUE){
    ## Get mean expression values of replicate samples from the same drug treatment in cells
    pert_cell_list <- split(inst$inst_id, inst$pert_cell_factor)
    chunk_list <- suppressWarnings(
        split(pert_cell_list, rep(seq_len(ceiling(length(pert_cell_list)/chunksize)), each=chunksize)))
    ## Creat an empty h5file
    if(file.exists(h5file)){
        if(isTRUE(overwrite)){
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