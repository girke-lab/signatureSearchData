#' Function processes the CEL files from each chip type (three for CMap2 data)
#' separately using the MAS5 normalization algorithm. The results will be written
#' to subdirectores (under a data parent directory) that are named after the chip
#' type names. To reduce the memory consumption of this step, the CEL files are
#' processed in user definable batch sizes. For details on the overall workflow,
#' please consult the vignette of this package. The normalization takes about 10
#' hours without parallelization. To save time, this process can be easily
#' accelerated on a computer cluster.
#' @title Normalize CEL Files
#' @param chiptype_list list, storing CEL files in each chiptype
#' @param batchsize number of CEL files to import and to normalize at once
#' @param rerun TRUE or FALSE, whether to run the function
#' @return Files storing normalized expression values in matrix, here one matrix file 
#' per chip type
#' @importFrom affy ReadAffy
#' @importFrom affy mas5
#' @importFrom affy exprs 
#' @examples 
#' # chiptype_list <- split(names(chiptype), as.character(chiptype))
#' normalizeCel(chiptype_list, rerun=FALSE)
#' @export
normalizeCel <- function(chiptype_list, batchsize=300, rerun=TRUE) {
    if(rerun) {
        for(i in names(chiptype_list)) {
            celfiles <- chiptype_list[[i]]
            cel_list <- suppressWarnings(
                split(celfiles, rep(seq_len(ceiling(length(celfiles)/batchsize)), each=batchsize)))
            chipdir <- paste0("./data/", i)
            if(dir.exists(chipdir)) unlink(chipdir, recursive=TRUE)
            dir.create(chipdir)
            mydir <- getwd()
            setwd(chipdir)
            ## Function to run MAS5 in batch
            f <- function(x, cel_list, batchsize) {
                dir.create(paste0("cellbatch_", x))
                mydata <- affy::ReadAffy(filenames=cel_list[[x]], celfile.path="../CEL")
                eset <- affy::mas5(mydata)
                # eset_pma <- affy::mas5calls(mydata) # Generates MAS 5.0 P/M/A calls.
                write.table(affy::exprs(eset), file=paste0("cellbatch_", x, "/mas5exprs.xls"), 
                            quote=FALSE, sep="\t", col.names = NA) 
                }
            lapply(names(cel_list), f, cel_list=cel_list, batchsize=batchsize)
            setwd(mydir)
        }
    } else {
        print("To execute function, set 'rerun=TRUE'")
    }
}
