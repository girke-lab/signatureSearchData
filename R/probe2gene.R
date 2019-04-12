#' probe set to gene level data
#' 
#' Transform expression values from probe set to gene level data. If genes are 
#' represented by several probe sets then their mean intensities are used. 
#' Probe sets not matching any gene are removed.
#' @param expr_df data.frame with expression values in probe set level
#' @param annot data.frame containing probe set id to gene entrez id 
#' mapping information
#' @return data.frame with expression values in gene level
#' @examples
#' expr_df <- data.frame(t1=runif(3), t2=runif(3), row.names=paste0("p", 1:3))
#' annot <- data.frame(ENTREZID=c("123","123","124"), 
#'                     SYMBOL=c("g1","g1","g2"),
#'                     row.names=paste0("p", 1:3))
#' gdf <- probe2gene(expr_df, annot)
#' @export
probe2gene <- function(expr_df, annot){
    myAnnot <- annot[as.character(annot[,"ENTREZID"]) != "NA",]
    expr_df <- expr_df[rownames(myAnnot),]
    idlist <- tapply(row.names(myAnnot), as.character(myAnnot$ENTREZID), c)
    gexpr_df <- t(sapply(names(idlist), function(x) colMeans(expr_df[idlist[[x]], ])))
    return(gexpr_df)
}