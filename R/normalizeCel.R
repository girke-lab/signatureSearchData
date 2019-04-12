#' Process each chip type separately running MAS5 with BiocParallel on cluster
#' @title normalize CEL files
#' @param chiptype_list list, storing CEL files in each chiptype
#' @param partition partition of a computer cluster 
#' @param rerun TRUE or FALSE, whether to run the function
#' @return files storing normalized matrix in each chip type
#' @importFrom affy ReadAffy
#' @importFrom affy mas5
#' @importFrom affy mas5calls
#' @importFrom Biobase assayDataElement
#' @import batchtools
#' @examples 
#' # chiptype_list <- split(names(chiptype), as.character(chiptype))
#' normalizeCel(chiptype_list, rerun=FALSE)
#' @export
normalizeCel <- function(chiptype_list, partition="girkelab", rerun=TRUE) {
    if(rerun==TRUE) {
        for(i in names(chiptype_list)) {
            saveRDS(chiptype_list[[i]], "./data/chiptype_tmp.rds")
            celfiles <- readRDS("./data/chiptype_tmp.rds")
            batchsize <- 100
            cel_list <- suppressWarnings(
                split(celfiles, rep(seq_len(ceiling(length(celfiles)/batchsize)), each=batchsize)))
            dir.create(paste0("./data/", i))
            batchtools_file <- system.file("extdata", "batchtools.conf.R", package="signatureSearchData")
            slurm_file <- system.file("extdata", "slurm.tmpl", package="signatureSearchData")
            file.copy(c(slurm_file, batchtools_file), paste0("./data/", i))
            file.rename(paste0("./data/", i,"/batchtools.conf.R"), paste0("./data/", i,"/.batchtools.conf.R"))
            mydir <- getwd()
            setwd(paste0("./data/", i))
            ## Function to run MAS5 on cluster with BiocParallel
            f <- function(x) {
                celfiles <- readRDS("../chiptype_tmp.rds")
                batchsize <- 100
                cel_list <- suppressWarnings(
                    split(celfiles, rep(seq_len(ceiling(length(celfiles)/batchsize)), each=batchsize)))
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
            if(dir.exists("mas5_reg")) 
                unlink("mas5_reg", recursive=TRUE)
            reg = makeRegistry(file.dir="mas5_reg", conf.file=".batchtools.conf.R") 
            ids = batchMap(f, x = names(cel_list))
            done <- submitJobs(ids, resources=list(walltime=600, ncpus=1, memory=10240, 
                                                   partition=partition), reg=reg) 
            if(!waitForJobs()) warning(paste0("Not all jobs are successfully completed under ",
                                              i, ", please check!")) 
            setwd(mydir)
            unlink("./data/chiptype_tmp.rds")
        }
    } else {
        print("To execute function, set 'rerun=TRUE'")
    }
}

#' Process each chip type separately running rma normalization
#' 
#' The rma method normalize all the CEL files in the same chip type as one batch,
#' it need to allocate large number of memory, 512Gb works
#' @title normalize CEL files with rma method
#' @param chiptype_list list, storing CEL files in each chiptype
#' @param rerun TRUE or FALSE, whether to run the function
#' @return files storing rma normalized matrix in each chip type
#' @importFrom affy ReadAffy
#' @importFrom affy rma
#' @importFrom affy exprs
#' @examples 
#' # chiptype_list <- split(names(chiptype), as.character(chiptype))
#' normalizeCel_rma(chiptype_list, rerun=FALSE)
#' @export
normalizeCel_rma <- function(chiptype_list, rerun = FALSE){
    mydir <- getwd()
    if(rerun==TRUE){
        for(i in names(chiptype_list)){
            if(! dir.exists(paste0("./data/",i)))
                dir.create(paste0("./data/",i))
            setwd(paste0("./data/",i))
            mydata <- affy::ReadAffy(filenames=chiptype_list[[i]], 
                                     celfile.path="../CEL")
            eset <- affy::rma(mydata)
            write.table(affy::exprs(eset), file="./rma_exprs.xls", 
                        quote=FALSE, sep="\t",col.names=NA)
            setwd("../../")
        }
    } else {
        print("To execute function, set 'rerun=TRUE'")
    }
    setwd(mydir)
}
