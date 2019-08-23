#' This function generates a list of treatment samples and control samples CEL ids
#' for each compound treatment from the experiment annotation of CMAP02 
#' instances. For more context of the CMAP02 database, please consult the vignette
#' of this package. 
#' @title Generate List of Treatment vs. Control Samples 
#' @param cmap data.frame containing experiment information of CMAP02 instances
#' @param myby "CMP" or "CMP_CELL", "CMP": by compound treatments in all cells;
#' "CMP_CELL": by compound treatments in individual cell
#' @return a list object containing treatment samples and control samples ids
#' for each compound treament. 
#' @examples 
#' path <- system.file("extdata", "cmap_instances_02.txt", package="signatureSearchData")
#' cmap_inst <- read.delim(path, check.names=FALSE) 
#' comp_list <- sampleList(cmap_inst, myby="CMP_CELL")
#' @export
sampleList <- function(cmap, myby) {
    ## Reconstruct CEL file names for control samples
    byCel_list <- split(cmap[, c("perturbation_scan_id", "vehicle_scan_id4")], 
                        as.character(cmap$perturbation_scan_id))
    sampleList1 <- function(x) {
        c <- unlist(strsplit(as.character(x[,2]), "\\."))
        if(max(nchar(c)>4)) {    
            c <- gsub("'", "", as.character(x[,2]))
        } else {
            c <- c[nchar(c)!=0]
            c <- paste0(gsub("\\..*", "", as.character(x[,1])), ".", c)
            c <- gsub("'", "", c) # Quite a few CEL file names have a "'" prepended 
        }
        t <- gsub("'", "", as.character(x[,1]))
        return(list(list(t=t, c=c)))
    }
    byCel_list <- sapply(byCel_list, sampleList1)
    
    ## Split by CMP only, or CMP and cell type
    if(myby=="CMP") {
        byCMP_list <- split(cmap[, c("cmap_name", "perturbation_scan_id")], 
                            as.character(cmap$cmap_name))
    }
    if(myby=="CMP_CELL") {
        cmap <- data.frame(cmp_cell=paste(cmap$cmap_name, cmap$cell2, sep="_"), cmap)
        byCMP_list <- split(cmap[, c("cmap_name", "perturbation_scan_id")], 
                            as.character(cmap$cmp_cell))
    }
    sampleList2 <- function(x, y=byCel_list) {
        s <- y[x[,2]]    
        l <- list(t=paste0(unique(as.character(unlist(sapply(seq_along(s), function(x) s[[x]]$t)))), ".CEL"),
                  c=paste0(unique(as.character(unlist(sapply(seq_along(s), function(x) s[[x]]$c)))), ".CEL"))
        return(list(l))
    }
    byCel_list <- sapply(names(byCMP_list), function(i) sampleList2(x=byCMP_list[[i]], y=byCel_list))
    return(byCel_list)
}