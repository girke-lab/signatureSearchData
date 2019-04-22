#' Combine batched normalization results from same chip type into a single data frame
#' @title combine batched normalization results
#' @param chiptype_dir chiptype directories
#' @param rerun TRUE or FALSE, whether to run the function
#' @return file storing normalization values of CELs for each chiptype
#' @importFrom utils read.delim
#' @importFrom utils write.table
#' @examples 
#' # chiptype_dir <- unique(readRDS("./data/chiptype.rds"))
#' combineResults(chiptype_dir, rerun=FALSE)
#' @export
combineResults <- function(chiptype_dir, rerun=TRUE) {
    if(isTRUE(rerun)) {
        for(j in seq_along(chiptype_dir)) {
            mydirs <- list.files(paste0("data/", chiptype_dir[j]), pattern="cellbatch_", full.names=TRUE)
            for(i in seq_along(mydirs)) {
                if(i==1) {
                    df1 <- read.delim(paste0(mydirs[i], "/", "mas5exprs.xls"), row.names=1, check.names=FALSE)
                } else {
                    tmpdf1 <- read.delim(paste0(mydirs[i], "/", "mas5exprs.xls"), row.names=1, check.names=FALSE)
                    df1 <-cbind(df1, tmpdf1)
                }
                cat("Processed", i, "of", length(mydirs), "\n")
            }
            write.table(df1, paste0("data/", chiptype_dir[j], "/all_mas5exprs.xls"), 
                        quote=FALSE, sep="\t", col.names = NA) 
            saveRDS(df1, paste0("data/", chiptype_dir[j], "/", "all_mas5exprs.rds")) # For fast loading
            cat("Generated", paste0("data/", chiptype_dir[j], "/all_mas5exprs.xls"), "\n")
        }
    } else {
        print("To execute function, set 'rerun=TRUE'")
    }
}