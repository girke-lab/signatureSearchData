#' Combine normalization results from same chip type in single data frame
#' @title combine normalization results
#' @param chiptype_dir chiptype directories
#' @param rerun TRUE or FALSE, whether to run the function
#' @return files storing normalization values for each chiptype
#' @examples 
#' # chiptype_dir <- unique(readRDS("./data/chiptype.rds"))
#' combineResults(chiptype_dir, rerun=FALSE)
combineResults <- function(chiptype_dir, rerun=TRUE) {
    if(isTRUE(rerun)) {
        for(j in seq_along(chiptype_dir)) {
            mydirs <- list.files(paste0("data/", chiptype_dir[j]), pattern="cellbatch_", full.names=TRUE)
            for(i in seq_along(mydirs)) {
                if(i==1) {
                    df1 <- read.delim(paste0(mydirs[i], "/", "mas5exprs.xls"), row.names=1, check.names=FALSE)
                    df2 <- read.delim(paste0(mydirs[i], "/", "mas5pma.xls"), row.names=1, check.names=FALSE)
                    df3 <- read.delim(paste0(mydirs[i], "/", "mas5pval.xls"), row.names=1, check.names=FALSE)
                } else {
                    tmpdf1 <- read.delim(paste0(mydirs[i], "/", "mas5exprs.xls"), row.names=1, check.names=FALSE)
                    df1 <-cbind(df1, tmpdf1)
                    tmpdf2 <- read.delim(paste0(mydirs[i], "/", "mas5pma.xls"), row.names=1, check.names=FALSE)
                    df2 <-cbind(df2, tmpdf2)
                    tmpdf3 <- read.delim(paste0(mydirs[i], "/", "mas5pval.xls"), row.names=1, check.names=FALSE)
                    df3 <-cbind(df3, tmpdf3)
                }
                cat("Processed", i, "of", length(mydirs), "\n")
            }
            write.table(df1, paste0("data/", chiptype_dir[j], "/all_mas5exprs.xls"), 
                        quote=FALSE, sep="\t", col.names = NA) 
            saveRDS(df1, paste0("data/", chiptype_dir[j], "/", "all_mas5exprs.rds")) # For fast loading
            cat("Generated", paste0("data/", chiptype_dir[j], "/all_mas5exprs.xls"), "\n")
            write.table(df2, paste0("data/", chiptype_dir[j], "/all_mas5pma.xls"), 
                        quote=FALSE, sep="\t", col.names = NA) 
            cat("Generated", paste0("data/", chiptype_dir[j], "/all_mas5pma.xls"), "\n")
            write.table(df3, paste0("data/", chiptype_dir[j], "/all_mas5pval.xls"), 
                        quote=FALSE, sep="\t", col.names = NA) 
            cat("Generated", paste0("data/", chiptype_dir[j], "/all_mas5pval.xls"), "\n")
        }
    } else {
        print("To execute function, set 'rerun=TRUE'")
    }
}