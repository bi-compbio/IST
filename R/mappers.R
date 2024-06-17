#' @title Orthology mappers
#' 
#' @description \code{orthIDdouble} 
#' identifies ortholog genes using biomart homology table.
#' This function is based on YS's \code{orthIDcon}
#' from the \code{SCISSORS} repository
#' 
#' @param org.from character of species name 
#' (must match species name used in Biomart)
#' @param org.to character of species name 
#' (must match species name used in Biomart)
#' @param ... further arguments passed to \code{biomaRt::useMart()}, anything
#' besides \code{biomart} and \code{dataset} 
#' 
#' @return data.table of id conversion table
#' 
#' @examples
#' \dontrun{df.mouse2human <- orthIDdouble("mmusculus", "hsapiens")}
#'  
#' @name orthology-mappers
#' @rdname orthology-mappers
#' 
#' @importFrom checkmate qassert
#' @importFrom biomaRt useMart getBM
#' @import data.table
#' @export
orthIDdouble <- function(org.from = NULL, org.to = NULL, ...) {
    checkmate::qassert(org.from, "S1")
    checkmate::qassert(org.to, "S1")
    
    ensembl <- biomaRt::useMart(
        biomart = "ensembl", 
        dataset = paste0(org.from, "_gene_ensembl"), 
        ...)
    
    homoAttribute <- c(
        "homolog_ensembl_gene", 
        "homolog_orthology_confidence", 
        "homolog_orthology_type")
    
    homoAttribute <- paste0(org.to, "_", homoAttribute)
    myfilter <- paste0("with_", org.to, "_homolog")
    idcon <- biomaRt::getBM(
        attributes = c('ensembl_gene_id', homoAttribute),
        filter = myfilter,
        values = TRUE,
        mart = ensembl)
    
    idcon <- data.table::data.table(idcon)
    names(idcon) <- c("gene", "ortholog", "confidence", "type")
    
    idcon <- idcon[type == "ortholog_one2one", ]
    # filter out m to m relationship
    idcon <- idcon[
        !gene %in% idcon[, .N, by = gene][N > 1, gene] & 
            !ortholog %in% idcon[, .N, by = ortholog][N > 1, ortholog],
        ]
    
    idcon
}

#' @description \code{orthIDsingle} build a trivial homolgoy 
#' mapping in the within-species case
#' 
#' @param org character of species name 
#' (must match species name used in Biomart)
#' 
#' @examples
#' \dontrun{df.mouse <- orthIDsingle("mmusculus")}
#' 
#' @name orthology-mappers
#' @rdname orthology-mappers
#' 
#' @importFrom data.table data.table
#' @importFrom biomaRt useMart getBM
#' @export
orthIDsingle <- function(org = NULL, ...) {
    ensembl <- biomaRt::useMart(
        "ensembl", 
        dataset = paste0(org, "_gene_ensembl"), 
        ...)
    
    # data frame with one column
    idcon <- biomaRt::getBM(mart = ensembl, attributes = "ensembl_gene_id")
    
    data.table::data.table(
        gene = idcon[, 1], 
        ortholog = idcon[, 1], 
        confidence = 1, 
        type = "same_gene"
    )
}

#' @description \code{orthIDcon} covers the general case, dispatching 
#' the adequate function and setting proxy settings
#' 
#' @param proxy character, proxy for the database connection; leave missing 
#' for not setting any proxy at all
#' 
#' @examples
#' \dontrun{
#' df.mouse <- orthIDcon(org.from = "mmusculus", org.to = "mmusculus")
#' df.mmu2hsa <- orthIDcon(org.from = "mmusculus", org.to = "hsapiens")
#' }
#' 
#' @name orthology-mappers
#' @rdname orthology-mappers
#' 
#' @importFrom data.table data.table
#' @importFrom biomaRt useMart getBM
#' @export
orthIDcon <- function(org.from = NULL, org.to = NULL, proxy, ...) {
    # set proxy if desired
    if (!missing(proxy)) {
        checkmate::qassert(proxy, "S1")
        Sys.setenv(http_proxy = proxy)
    }
    
    # same org.from and org.to
    if (org.from == org.to) {
        orthIDsingle(org = org.to, ...)
    } else {
        orthIDdouble(org.from = org.from, org.to = org.to, ...)
    }
}