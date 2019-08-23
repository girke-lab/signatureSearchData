#' This function runs DEG analysis with Limma at user defined FDR and LFC
#' cutoff by providing CMAP02 normalized expression values and annotation of
#' treatment v.s. control samples. For more context of the CMAP02 database, 
#' please consult the vignette of this package. 
#' @title DEG Analysis with Limma
#' @param df data.frame containing normalized intensity values of CMAP02 samples
#' @param comp_list list of CEL file ids of treatment and control samples for 
#' each compound treatment. The list for CMAP02 data is generated from 
#' the \code{\link{sampleList}} function.
#' @param fdr cutoff of false discovery rate (FDR) for defining DEGs
#' @param foldchange cutoff of log2 fold change (LFC) for defining DEGs
#' @param verbose TRUE or FALSE
#' @return list containing DEGs and log2FC matrix
#' @importFrom Biobase ExpressionSet
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom stats model.matrix
#' @examples 
#' path <- system.file("extdata", "cmap_instances_02.txt", package="signatureSearchData")
#' cmap_inst <- read.delim(path, check.names=FALSE) 
#' comp_list <- sampleList(cmap_inst, myby="CMP_CELL")
#' df <- as.data.frame(matrix(runif(70), ncol=7))
#' colnames(df) <- unlist(comp_list[1])
#' degList <- runLimma(df, comp_list[1])  
#' @export
runLimma <- function(df, comp_list, fdr=0.05, foldchange=1, verbose=TRUE) {
    ## Generate result container
    deg <- matrix(0, nrow=nrow(df), ncol=length(comp_list), 
                  dimnames = list(rownames(df), names(comp_list)))
    colnames(deg) <- names(comp_list)
    #deg_list <- NULL # only used if affyid not NULL
    logfc_ma <- deg
    pvalue_ma <- deg
    ## Run limma
    for(i in seq_along(comp_list)) {
        sample_set <- unlist(comp_list[[i]])
        repcounts <- sapply(comp_list[[i]], length)
        repcounts <- paste("rep count:", paste(paste0(names(repcounts), repcounts), 
                                               collapse="_"))
        ## The following if statement will skip sample sets with less than 3 CEL files 
        ## (control and treatment) or those containing non-existing CEL files. 
        ## For tracking NAs will be injected in the corresponding columns of the
        ## result matrix.
        if((length(sample_set)<3) | any(!sample_set %in% colnames(df))) {
            deg[, i] <- NA
            if(verbose==TRUE) {
                cat("Sample", i, "of", paste0(length(comp_list), ":"), 
                    paste0("not enough usable replicates (", repcounts, ")."), "\n")
            }
        } else {
            dfsub <- df[, sample_set]
            # eset <- new("ExpressionSet", exprs = as.matrix(dfsub), annotation="hgu133a")
            eset <- Biobase::ExpressionSet(assayData=as.matrix(dfsub),  annotation="hgu133a")
            repno <- rep(seq_along(comp_list[[i]]), sapply(comp_list[[i]], length))
            design <- model.matrix(~ -1+factor(repno))
            colnames(design) <- names(comp_list[[i]])
            # Fit a linear model for each gene based on the given series of arrays
            fit <- limma::lmFit(eset, design) 
            contrast.matrix <- limma::makeContrasts(contrasts="t-c", levels=design)
            fit2 <- limma::contrasts.fit(fit, contrast.matrix)
            # Computes moderated t-statistics and log-odds of differential expression 
            # by empirical Bayes shrinkage of the standard errors towards a common value.
            fit2 <- limma::eBayes(fit2)
            limmaDF <- limma::topTable(fit2, coef=1, adjust="fdr", sort.by="none", number=Inf)
            # if(!is.null(affyid)) {
            #     tmp_list <- list(limmaDF[affyid,])
            #     names(tmp_list) <- names(comp_list[i])
            #     deg_list <- c(deg_list, tmp_list)
            # } 
            pval <- limmaDF$adj.P.Val <= fdr # FDR 1%
            # Fold change 2
            fold <- (limmaDF$logFC >= foldchange | limmaDF$logFC <= -foldchange) 
            affyids <- rownames(limmaDF[pval & fold,])
            deg[affyids, i] <- 1
            pvalue_ma[ , i] <- limmaDF[, "adj.P.Val"]
            logfc_ma[ , i] <- limmaDF[,"logFC"]
            if(verbose==TRUE) {
                cat("Sample", i, "of", paste0(length(comp_list), ":"), "identified", 
                    length(affyids), paste0("DEGs (", repcounts, ")."), "\n")
            }
        }
    }
    # if(!is.null(affyid)) {
    #     return(deg_list)
    #} 
    col_index <- !is.na(colSums(deg)) # Remove columns where DEG analysis was not possible
    return(list(DEG=deg[, col_index], logFC=logfc_ma[, col_index], FDR=pvalue_ma[, col_index]))
}
