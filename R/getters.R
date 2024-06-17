#' @title Get class slots and processed tables
#' 
#' @description \code{getTabDecisions()} returns the decisions of the specified
#' type of model \code{mod}
#' 
#' @param ist.results object holding the results (must be fitted)
#' @param mod character, type of model (oc, bin)
#' 
#' @return data.table with the desired data
#' 
#' @name getters
#' @rdname getters
#' 
#' @importFrom methods slot
#' @importFrom checkmate assertClass assertSubset qassert
#' @export
getTabDecisions <- function(ist.results, mod = "bin", ...) {
    checkmate::assertClass(ist.results, "ist.results")
    checkmate::qassert(mod, "S1")
    checkmate::assertSubset(mod, c("oc", "bin"))
    
    methods::slot(ist.results, paste0("tab.decisions.", mod))
}

#' @description \code{getTabMetrics()} returns the performance
#' metrics of the specified type of model \code{mod}
#' 
#' @rdname getters
#' 
#' @importFrom methods slot
#' @importFrom checkmate assertClass assertSubset qassert
#' @export
getTabMetrics <- function(ist.results, mod = "bin", ...) {
    checkmate::assertClass(ist.results, "ist.results")
    checkmate::qassert(mod, "S1")
    checkmate::assertSubset(mod, c("oc", "bin"))
    
    methods::slot(ist.results, paste0("tab.metrics.", mod))
}

#' @description \code{getTabPls()} returns the coefficients and 
#' scale parameters in the PLS models in BIN
#' 
#' @rdname getters
#' 
#' @importFrom checkmate assertClass
#' @export
getTabPls <- function(ist.results, ...) {
    checkmate::assertClass(ist.results, "ist.results")
    
    ist.results@tab.pls.bin
}

#' @description \code{getTabSignatures()} returns the 
#' signatures in tabular format from an \code{ist.results} object
#' 
#' @rdname getters
#' 
#' @importFrom checkmate assertClass
#' @export
getTabSignatures <- function(ist.results, ...) {
    checkmate::assertClass(ist.results, "ist.results")
    
    ist.results@tab.signatures
}

#' @description \code{getTabDelta()} returns the delta 
#' tables of an ist.results object 
#' at the median (reference data), gene and pathway level
#' (using the \code{level} argument)
#' 
#' @rdname getters
#' 
#' @importFrom methods slot
#' @importFrom checkmate assertClass assertSubset qassert
#' @export
getTabDelta <- function(ist.results, level = "pathway", ...) {
    checkmate::assertClass(ist.results, "ist.results")
    checkmate::qassert(level, "S1")
    checkmate::assertSubset(level, c("median", "gene", "pathway"))
    
    methods::slot(ist.results, paste0("tab.delta.", level))
}

#' @description \code{getMetaSig()} returns the metadata table for the 
#' signatures
#' 
#' @rdname getters
#' 
#' @export
getMetaSig <- function(object, ...) {

    if (inherits(object, "ist.results")) {
        obj <- object@ist.signatures
    } else if (inherits(object, "ist.signatures")) {
        obj <- object
    } else {
        stop("object is not of class ist.results or ist.signatures")
    }
    
    obj@tab.meta
}

#' @description \code{getMetaPathways()} returns the metadata table for the 
#' pathways
#' 
#' @rdname getters
#' 
#' @export
getMetaPathways <- function(object, mod = "bin", ...) {
    checkmate::qassert(mod, "S1")
    checkmate::assertSubset(mod, c("oc", "bin"))
    
    if (inherits(object, "ist.results")) {
        obj <- object@ist.pathways
    } else if (inherits(object, "ist.signatures")) {
        obj <- object
    } else {
        stop("object is not of class ist.results or ist.pathways")
    }
    
    arg.mod <- mod
    obj@pathways.meta[mod == arg.mod, ]
}

#' @description \code{getPathways()} returns the pathway is in an 
#' \code{ist.pathways} or an \code{ist.results} object
#' 
#' @param object \code{ist.pathways}, \code{ist.results} or 
#' \code{ist.signatures}
#' @param mod character, type of model (oc, bin)
#' 
#' @return getPathways, getSignatures: character vector with names
#' 
#' @name getters
#' @rdname getters
#' 
#' @importFrom methods slot
#' @importFrom checkmate assertClass assertSubset qassert
#' @export
getPathways <- function(object, mod = "bin", ...) {

    checkmate::qassert(mod, "S1")
    checkmate::assertSubset(mod, c("oc", "bin"))
    
    if (inherits(object, "ist.results")) {
        obj <- object@ist.pathways
    } else if (inherits(object, "ist.pathways")) {
        obj <- object
    } else {
        stop("object is not of class ist.results or ist.pathways")
    }
    
    slot.nm <- paste0("list.", mod)
    names(methods::slot(obj, slot.nm))
}

#' @description \code{getSignatures()} returns the signature ids in an 
#' \code{ist.pathways} or an \code{ist.signatures} object
#' 
#' @name getters
#' @rdname getters
#' 
#' @export
getSignatures <- function(object, ...) {
    
    if (inherits(object, "ist.results")) {
        obj <- object@ist.signatures
    } else if (inherits(object, "ist.signatures")) {
        obj <- object
    } else {
        stop("object is not of class ist.results or ist.signatures")
    }
    
    names(obj@list.data)
}

#' @description \code{getGroups()} returns the sample groups in an 
#' \code{ist.results} object, whereas \code{getGroupLevels()}
#' returns their levels
#' 
#' @name getters
#' @rdname getters
#' 
#' @importFrom checkmate assertClass
#' @export
getGroups <- function(object, ...) {
    
    checkmate::assertClass(object, "ist.results")
    
    object@group
}

#' @importFrom checkmate assertClass 
#' @export
getGroupLevels <- function(object, ...) {
    
    checkmate::assertClass(object, "ist.results")
    
    levels(object@group)
}

#' @description \code{getDataMatrix()} returns the numeric matrix with 
#' the readouts that are used to train the models in an \code{ist.pathways}
#' object, or an \code{ist.results} containing it, and 
#' \code{getDataResponse()} gives the \code{y} vector
#' 
#' @name getters
#' @rdname getters
#' 
#' @importFrom checkmate assertClass
#' @export
getDataMatrix <- function(object, ...) {
    
    if (inherits(object, "ist.results")) {
        obj <- object@ist.pathways
    } else if (inherits(object, "ist.pathways")) {
        obj <- object
    } else {
        stop("object is not of class ist.results or ist.pathways")
    }
    
    obj@X
}

#' @name getters
#' @rdname getters
#' 
#' @importFrom checkmate assertClass
#' @export
getDataResponse <- function(object, ...) {
    
    if (inherits(object, "ist.results")) {
        obj <- object@ist.pathways
    } else if (inherits(object, "ist.pathways")) {
        obj <- object
    } else {
        stop("object is not of class ist.results or ist.pathways")
    }
    
    obj@y
}

#' @description \code{getSampleIds()} returns the sample ids, whereas
#' \code{getGeneIds()} returns the gene identifiers from the readouts 
#' that were used to train the model 
#' (\code{ist.pathways} or \code{ist.results})
#' 
#' @name getters
#' @rdname getters
#' 
#' @importFrom checkmate assertClass
#' @export
getSampleIds <- function(object, ...) {
    rownames(getDataMatrix(object))
}

#' @name getters
#' @rdname getters
#' 
#' @importFrom checkmate assertClass
#' @export
getGeneIds <- function(object, ...) {
    colnames(getDataMatrix(object))
}

#' @title Extract pls weights from model in a data.table
#' 
#' @param object an ist.mod.bin object
#' 
#' @importFrom stats coef
#' @import data.table
get.plstab.bin <- function(object) {
    # do not fail in empty models
    if (is.null(object)) return(NULL)
    
    pls.coef <- stats::coef(object@pls.mod, ncomp = object@Ncomp)[, , 1]
    pls.scale <- object@pls.mod$scale
    
    df.pls <- data.frame(
        ortholog = names(pls.coef), 
        Coef = pls.coef, 
        Scale = pls.scale
    )
    
    data.table::as.data.table(df.pls)
}

#' @title Extract performance metrics from ist.mod.oc object
#' 
#' @param object an ist.mod.oc object
#' 
#' @import data.table
get.perf.oc <- function(object) {
    if (is.null(object)) return(NULL)
    
    ans <- data.table::as.data.table(object@svm.tuned$performances)
    data.table::setnames(
        ans, 
        old = c("error", "dispersion"), 
        new = c("Accuracy", "AccuracySD"))
    ans
}

#' @title Extract performance metrics from ist.mod.bin object
#' 
#' @param object an ist.mod.bin object
#' 
#' @import data.table
get.perf.bin <- function(object) {
    if (is.null(object)) return(NULL)
    
    data.table::as.data.table(object@train$results)
}

#' @title Extract table to plot boxplots
#' 
#' @param object an ist.results
#' @param quote.sub an expression enclosed by quote(), to 
#' subset signatures (sig.id) and pathways (pathway). 
#' If left to TRUE, it selects all signatures and pathways
#' 
#' @import data.table
get.tab.boxplot <- function(object, quote.sub = quote(TRUE)) {

    tab.boxplot <- getTabDecisions(object)[eval(quote.sub)]
    getMetaSig(object)[tab.boxplot, on = "sig.id"]
}


#' @title Extract table to plot gene heatmaps
#'
#' @param object an ist.results
#' @param id.path character, pathway to represent
#' @param sig.ids character, signature ids to include (default: all)
#' @param long.only logical, whether to return the long table 
#' (ggplot) or the wide data with row/column annotations (pheatmap)
#' @param vars.meta.sig character vector, column names from the signature 
#' metadata to include in the graphical annotations. Mandatory if 
#' \code{long.only = FALSE}
#' @param max.genes numeric, maximum number of genes to display in genemap
#' (defaults to all). Genes will be prioritised using their sum of squares
#'
#' @import data.table
get.tab.genemap <- function(
    object, id.path, sig.ids = getSignatures(x), 
    long.only = TRUE, vars.meta.sig, max.genes = Inf) {

    # get gene contribution tables, subset to pathway and signature
    df.hm <- getTabDelta(object, "gene")
    df.hm <- df.hm[pathway == id.path & sig.id %in% sig.ids]
    
    if (nrow(df.hm) == 0) {
        message(
            "No overlap between genes in pathway ", id.path, 
            " and genes in the signatures in sig.id. ", 
            "Returning NULL genemap")
        return(NULL)
    } 
    
    # if a gene mapping is provided, change the labels
    gene2label <- object@vec.gene2label
    if (length(gene2label) > 0) {
        message("Using custom gene labels")
        df.hm[, ortholog := gene2label[as.character(ortholog)]]
    }
    
    # get pls row annotations
    df.pls <- getTabPls(object)[pathway == id.path]
    df.pls[, Weight := Coef/Scale]    
    if (length(gene2label) > 0) {
        df.pls[, ortholog := gene2label[as.character(ortholog)]]
    }

    # wide format matrix
    df.wide.hm <- data.table::dcast(
        df.hm[!is.na(sig.id)], sig.id ~ ortholog,
        value.var = "delta.percent", fill = 0)

    mat.hm <- as.matrix(df.wide.hm[, -c("sig.id")])
    rownames(mat.hm) <- df.wide.hm$sig.id
    
    # leave low-weight genes out of heatmap if max.genes imposes so
    #   add a sort() so that the original column order is kept (just in case)
    mat.sqsum <- colSums(mat.hm**2)
    mat.colkeep <- sort(tail(order(mat.sqsum), max.genes))
    mat.hm <- mat.hm[, mat.colkeep, drop = FALSE]

    # columns
    col.hm <- as.data.frame(df.pls)
    rownames(col.hm) <- df.pls$ortholog
    # only take the genes in the heatmap
    col.hm <- col.hm[
        colnames(mat.hm), c("Weight", "Coef", "Scale"), drop = FALSE]

    # rows
    if (missing(vars.meta.sig)) vars.meta.sig <- NULL
    if (!is.null(vars.meta.sig)) {
        row.hm <- as.data.frame(getMetaSig(object))
        rownames(row.hm) <- row.hm$sig.id
        row.hm <- row.hm[rownames(mat.hm), vars.meta.sig, drop = FALSE]
    } else {
        row.hm <- NULL
    }
    
    # subset long format too
    df.hm <- df.hm[ortholog %in% rownames(col.hm)]
    
    # return statements
    if (long.only) return(df.hm)
    
    list(
        row = row.hm, col = col.hm, 
        data.wide = mat.hm, data.long = df.hm)
}

#' @title Extract table to plot pathway heatmaps
#'
#' @param object an ist.results
#' @param id.path character, pathways to represent (default: all)
#' @param sig.ids character, signature ids to include (default: all)
#' @param long.only logical, whether to return the long table 
#' (ggplot) or the wide data with row/column annotations (pheatmap)
#' @param vars.meta.sig character vector, column names from the signature 
#' metadata to include in the graphical annotations. Mandatory if 
#' \code{long.only = FALSE}
#'
#' @import data.table
get.tab.pathwaymap <- function(
    object, id.path = getPathways(object), sig.ids = getSignatures(object), 
    long.only = TRUE, vars.meta.sig, vars.meta.path) {
    
    df.mhm <- getTabDelta(object, "pathway")
    df.mhm <- df.mhm[pathway %in% id.path & sig.id %in% sig.ids]
    
    if (nrow(df.mhm) == 0) {
        message(
            "No overlap between genes in pathways in id.path ", 
            " and genes in the signatures in sig.id. ", 
            "Returning NULL pathwaymap")
        return(NULL)
    } 
    
    if (long.only) return(df.mhm)
    
    # matrix
    df.wide.mhm <- data.table::dcast(
        df.mhm, 
        sig.id ~ pathway, 
        value.var = "total.delta.percent", 
        fill = 0)
    mat.mhm <- as.matrix(df.wide.mhm[, !"sig.id"])
    rownames(mat.mhm) <- df.wide.mhm[, sig.id]
    
    # columns: used-defined annotation
    if (missing(vars.meta.path)) vars.meta.path <- NULL
    if (!is.null(vars.meta.path)) {
        col.mhm <- as.data.frame(getMetaPathways(object, "bin"))
        rownames(col.mhm) <- col.mhm$path.id
        col.mhm <- col.mhm[colnames(mat.mhm), vars.meta.path, drop = FALSE]
    } else {
        col.mhm <- NULL
    }
    
    # rows
    if (missing(vars.meta.sig)) vars.meta.sig <- NULL
    if (!is.null(vars.meta.sig)) {
        row.mhm <- as.data.frame(getMetaSig(object))
        rownames(row.mhm) <- row.mhm$sig.id
        row.mhm <- row.mhm[rownames(mat.mhm), vars.meta.sig, drop = FALSE]
    } else {
        row.mhm <- NULL
    }

    list(
        row = row.mhm, col = col.mhm,
        data.wide = mat.mhm, data.long = df.mhm)
}