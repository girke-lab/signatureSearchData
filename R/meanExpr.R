#' Take mean expression values of drug treatment samples
#' 
#' Take mean expression values of multiple drug treatment samples at different 
#' concentration and duration in the same cell as expression value of that drug
#' treatment in the cell.
#' @param expr_df data.frame with expression values in CEL samples
#' @param cmap_inst data.frame defining drug treatment of CEL samples
#' @return data.frame with mean expression values of drug treatment in a cell
#' @examples
#' path <- system.file("extdata", "cmap_instances_02.txt", package="signatureSearchData")
#' cmap_inst <- read.delim(path, check.names=FALSE) 
#' expr_df <- as.data.frame(matrix(runif(30), ncol=5))
#' colnames(expr_df) <- paste0(cmap_inst$perturbation_scan_id[1:5],".CEL")
#' mexpr <- meanExpr(expr_df, cmap_inst[1:5,])
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