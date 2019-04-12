#' Combine normalization results into a single data frame from difference
#' chip types
#' @title combine normalization results from different chip types
#' @param chiptype_dir chiptype directories
#' @param norm_method normalization method, one of "MAS5" or "rma".
#' @return data frame storing normalization values combined from chiptypes
#' @examples 
#' # chiptype_dir <- unique(readRDS("./data/chiptype.rds"))
#' combineNormRes(chiptype_dir, norm_method="not run")
#' @export
combineNormRes <- function(chiptype_dir, norm_method) {
    if(norm_method=="MAS5"){
        df1 <- readRDS(paste0("data/", chiptype_dir[1], "/", "all_mas5exprs.rds"))
        df2 <- readRDS(paste0("data/", chiptype_dir[2], "/", "all_mas5exprs.rds"))
        df3 <- readRDS(paste0("data/", chiptype_dir[3], "/", "all_mas5exprs.rds"))
    } else if (norm_method=="rma"){
        df1 <- read.delim(paste0("./data/", chiptype_dir[1], "/rma_exprs.xls"), 
                          sep="\t", header=TRUE, row.names=1, check.names=FALSE)
        df2 <- read.delim(paste0("./data/", chiptype_dir[2], "/rma_exprs.xls"), 
                          sep="\t", header=TRUE, row.names=1, check.names=FALSE)
        df3 <- read.delim(paste0("./data/", chiptype_dir[3], "/rma_exprs.xls"), 
                          sep="\t", header=TRUE, row.names=1, check.names=FALSE)
    } else {
        message("Please set norm_method as one of 'MAS5' or 'rma'!")
        return(NULL)
    }
    affyid <- rownames(df1)[rownames(df1) %in% rownames(df2)]
    affyid <- affyid[affyid %in% rownames(df3)]
    normdf <- cbind(df1[affyid,], df2[affyid,], df3[affyid,])
    return(normdf)
}
