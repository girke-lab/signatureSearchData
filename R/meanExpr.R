#' Calculate Mean Value for Replicated Samples
#' 
#' Function averages the normalized gene expression values for different
#' concentrations, treatment times and replicates of compounds and cell types.
#' For more context on this step, please consult the corresponding section in the
#' vignette of this package.
#' @param expr_df \code{data.frame} containing normalized expression values
#' @param cmap_inst \code{data.frame} defining experimental conditions
#' @return \code{data.frame} with mean expression values of replicates
#' @examples
#' path <- system.file("extdata", "cmap_instances_02.txt", package="signatureSearchData")
#' cmap_inst <- read.delim(path, check.names=FALSE) 
#' expr_df <- as.data.frame(matrix(runif(30), ncol=5))
#' colnames(expr_df) <- paste0(cmap_inst$perturbation_scan_id[seq_len(5)],".CEL")
#' mexpr <- meanExpr(expr_df, cmap_inst[seq_len(5),])
#' @export
meanExpr <- function(expr_df, cmap_inst){
    cmap_inst <- data.frame(cmap_inst, 
            drug_cell = paste(cmap_inst$cmap_name, cmap_inst$cell2, sep ="_"))
    celid <- paste0(cmap_inst$perturbation_scan_id,".CEL")
    celid2 <- gsub("'", "", celid) # get rid of "'" in celid
    drug_celid_list <- tapply(celid2, cmap_inst$drug_cell, c) 
    drug_cell_expr <- sapply(names(drug_celid_list), function(x) 
        rowMeans(as.data.frame(expr_df[ ,drug_celid_list[[x]]])))
    return(drug_cell_expr)
}