#' Process each chip type separately running MAS5 with BiocParallel on cluster
#' @title normalize CEL files
#' @param chiptype_list list, storing CEL files in each chiptype
#' @param rerun TRUE or FALSE, whether to run the function
#' @return files storing normalized matrix in each chip type
#' @examples 
#' # chiptype_list <- split(names(chiptype), as.character(chiptype))
#' normalizeCel(chiptype_list, rerun=FALSE)
normalizeCel <- function(chiptype_list, rerun=TRUE) {
    if(rerun==TRUE) {
        for(i in names(chiptype_list)) {
            saveRDS(chiptype_list[[i]], "./data/chiptype_tmp.rds")
            celfiles <- readRDS("./data/chiptype_tmp.rds")
            batchsize <- 100
            cel_list <- suppressWarnings(
                split(celfiles, rep(seq_len(ceiling(length(celfiles)/batchsize)), each=batchsize)))
            dir.create(paste0("./data/", i))
            batchjobs_file <- system.file("extdata", ".BatchJobs.R", package="signatureSearchData")
            torque_file <- system.file("extdata", "torque.tmpl", package="signatureSearchData")
            file.copy(c(torque_file, batchjobs_file), paste0("./data/", i))
            mydir <- getwd()
            setwd(paste0("./data/", i))
            ## Function to run MAS5 on cluster with BiocParallel
            f <- function(x) {
                celfiles <- readRDS("../chiptype_tmp.rds")
                batchsize <- 100
                cel_list <- suppressWarnings(
                    split(celfiles, rep(1:(ceiling(length(celfiles)/batchsize)), each=batchsize)))
                dir.create(paste0("cellbatch_", x))
                mydata <- affy::ReadAffy(filenames=cel_list[[x]], celfile.path="../CEL")
                eset <- affy::mas5(mydata)
                eset_pma <- affy::mas5calls(mydata) # Generates MAS 5.0 P/M/A calls.
                write.table(affy::exprs(eset), file=paste0("cellbatch_", x, "/mas5exprs.xls"), 
                            quote=FALSE, sep="\t", col.names = NA) 
                write.table(affy::exprs(eset_pma), file=paste0("cellbatch_", x, "/mas5pma.xls"), 
                            quote=FALSE, sep="\t", col.names = NA) 
                write.table(Biobase::assayDataElement(eset_pma, "se.exprs"), 
                            file=paste0("cellbatch_", x, "/mas5pval.xls"), 
                            quote=FALSE, sep="\t", col.names = NA) 
            }
            funs <- makeClusterFunctionsTorque("torque.tmpl")
            param <- BatchJobsParam(
                length(cel_list), 
                resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="12gb"), 
                cluster.functions=funs)
            register(param)
            resultlist <- bplapply(names(cel_list), f)
            setwd(mydir)
            unlink("./data/chiptype_tmp.rds")
        }
    } else {
        print("To execute function, set 'rerun=TRUE'")
    }
}