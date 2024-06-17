#' @title Overlay fold changes
#' 
#' @description \code{overlay.fc()} adds the fold changes of the 
#' matching genes (by name) to the columns (each column is a gene) of 
#' the input data matrix
#' 
#' @param mat matrix, samples in rows and features in columns, 
#' must have row and colnames
#' @param v.fc named vector with log fold changes to apply. 
#' These are added directly to the gene readouts.
#' Zeroes and missing genes imply no change.
#' 
#' @return A transformed matrix with the dimensions of \code{mat}
#' 
#' @importFrom stats setNames
overlay.fc <- function(mat, v.fc) {
    gene.mat <- colnames(mat)
    gene.intersect <- intersect(gene.mat, names(v.fc))
    n.intersect <- length(gene.intersect)
    
    # 0s in the unmapped genes, FC otherwise
    gene.fcs <- stats::setNames(numeric(length(gene.mat)), gene.mat)
    gene.fcs[gene.intersect] <- v.fc[gene.intersect]
    
    sweep(mat, 2, gene.fcs, "+")
}

#' @title Predict methods
#' 
#' @description \code{predict.ist.translator()} in-silico-treats 
#' data from org.to
#' using a previously fitted \code{ist.translator}. Currently just 
#' adds the fold changes to the matching genes.
#' 
#' @param object a fitted \code{ist.translator} for 
#' \code{predict.ist.translator()}, 
#' \code{ist.discriminator} for \code{predict.ist.discriminator()}, 
#' \code{ist.mod.oc} for \code{predict.ist.oc()}, and 
#' \code{ist.mod.bin} for \code{predict.ist.bin()}
#' @param newdata matrix, samples in rows and features in columns, 
#' must have names
#' @param type character, type of prediction (\code{"decision"}, 
#' \code{"scores"})
#' @param ... \code{predict.ist.translator()} ignores it. 
#' \code{predict.ist.discriminator()} passes it to 
#' \code{predict.ist.oc()} (which ignores it) and \code{predict.ist.bin()}
#' (which ignores it unless \code{type == "scores"}, passing 
#' \code{...} to \code{predict()} for \code{pls} objects)
#' 
#' @return \code{predict.ist.translator()} returns a matrix 
#' with the transformed signatures. 
#' @return \code{predict.ist.signatures()} returns a list of matrices 
#' with the transformed signatures, one per signature. 
#' \code{predict.ist.oc()} and \code{predict.ist.bin(type = "decision")} 
#' return a data.frame with the decision values.
#' \code{predict.ist.bin(type = "scores")} returns a data.frame with the 
#' principal components as columns.
#' \code{predict.ist.discriminator()} wraps around both models, 
#' returning a data.frame where \code{flavour} and \code{sample} are colnames
#' 
#' @name predict-ist
#' @rdname predict-ist
#' 
#' @importFrom checkmate assertMatrix assertClass
predict.ist.translator <- function(object, newdata, ...) {
    
    checkmate::assertMatrix(
        newdata, mode = "numeric", any.missing = FALSE, 
        row.names = "unique", col.names = "unique")
    
    checkmate::assertClass(object, "ist.translator")
    if (nrow(object@mapping.data.table) == 0) 
        stop("The ist.translator must be fitted first, use fit()")
    
    v.fc <- object@mapping.data
    if (length(v.fc) == 0) 
        warning("No mapping data available, did you fit() the translator?")
    
    overlay.fc(mat = newdata, v.fc = v.fc)
}

#' @description \code{predict.ist.signatures()} in-silico-treats 
#' data from org.to
#' using a previously fitted \code{predict.ist.signatures}. 
#' 
#' @rdname predict-ist
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom checkmate assertMatrix assertClass
predict.ist.signatures <- function(object, newdata, ...) {
    
    checkmate::assertMatrix(
        newdata, mode = "numeric", any.missing = FALSE, 
        row.names = "unique", col.names = "unique")
    
    checkmate::assertClass(object, "ist.signatures")
    if (length(object@list.data) == 0) 
        stop("The ist.signatures must be fitted first, use fit()")
    
    
    message("In silico treating data")
    BiocParallel::bplapply(
        object@list.data, 
        overlay.fc, 
        mat = newdata
    )
}

#' @description \code{predict.ist.oc()} applies a 
#' fitted \code{ist.mod.oc} 
#' to new data, matching colnames.
#' 
#' @rdname predict-ist
#' 
#' @importFrom checkmate assertMatrix assertSubset
#' @importFrom stats predict
predict.ist.oc <- function(object, newdata, ...) {
    
    # if (missing(newdata)) 
    #     stop("The 'newdata' argument cannot be missing")
    checkmate::assertMatrix(
        newdata, mode = "numeric", any.missing = FALSE, 
        row.names = "unique", col.names = "unique")
    
    # check that we input all the required features
    checkmate::assertSubset(object@genes.top, colnames(newdata))
    
    newx <- newdata[, object@genes.top, drop = FALSE]
    
    oc.pred <- stats::predict(object@svm.best, newx, decision.values = TRUE)
    
    data.frame(
        decision.value = attr(oc.pred, "decision.values")[, 1], 
        decision = oc.pred, 
        sample = rownames(newx), 
        sample.intrain = rownames(newx) %in% object@id.mod)

}

#' @description \code{predict.ist.bin()} applies a fitted 
#' \code{ist.mod.bin} 
#' to new data, matching colnames.
#' 
#' @rdname predict-ist
#' 
#' @importFrom checkmate assertMatrix assertSubset
#' @importFrom stats predict
predict.ist.bin <- function(object, newdata, type = "decision", ...) {
    
    checkmate::assertMatrix(
        newdata, mode = "numeric", any.missing = FALSE, 
        row.names = "unique", col.names = "unique")
    
    checkmate::qassert(type, "S1")
    checkmate::assertSubset(type, c("decision", "scores"))
    
    # check that we input all the required features
    checkmate::assertSubset(object@genes.mod, colnames(newdata))
    
    newx <- newdata[, object@genes.mod, drop = FALSE]
    
    # by default, return the response as decision values
    if (type == "decision") {
        bin.pred <- as.vector(
            stats::predict(
                object@pls.mod, 
                newdata = newx, 
                ncomp = object@Ncomp,
                type = "response")
        )
        
        data.frame(
            decision.value = as.vector(bin.pred), 
            decision = as.vector(bin.pred) > 0, 
            sample = rownames(newx), 
            sample.intrain = rownames(newx) %in% object@id.mod)
    } else if (type == "scores") {
        # otherwise, let type be user-defined
        # here we use the user-defined number of components
        # (we don't force it through ncomp)
        # as it's usually for plotting
        # 
        # Added a drop() because predict returns an array of shape (X,1,Y)
        # I guess because the response is univariate
        mat.pred <- drop(
            stats::predict(
                object@pls.mod, 
                newdata = newx, 
                type = "scores", ...)
        )
        colnames(mat.pred) <- paste0("Comp", seq_len(ncol(mat.pred)))
        
        data.frame(
            sample = rownames(mat.pred),
            as.data.frame(mat.pred)
        )
    }
}

#' @description \code{predict.ist.discriminator()} applies all the fitted 
#' discriminator models within the object to 
#' new data (in-silico-treated or not)
#' 
#' @param ncomp numeric vector, only for binary classifiers (isd, isa) 
#' and \code{type = "scores"}, which components to return? 
#' Passed to \code{predict} on the \code{mvr} object
#' 
#' @rdname predict-ist
#' 
#' @importFrom methods slot
#' @importFrom plyr ldply
#' @importFrom checkmate assertMatrix assertClass qassert
#' @importFrom stats setNames
predict.ist.discriminator <- function(
    object, 
    newdata, 
    type = "decision",
    ncomp = c(1, 2),
    ...) {
    checkmate::assertMatrix(
        newdata, mode = "numeric", any.missing = FALSE, 
        row.names = "unique", col.names = "unique")
    
    checkmate::assertClass(object, "ist.discriminator")
    
    checkmate::qassert(type, "S1")
    checkmate::assertSubset(type, c("decision", "scores"))
    
    checkmate::qassert(ncomp, "N+")
    
    if (type == "decision") {
        ist.flavours <- c("iss", "isd", "isa")
    } else if (type == "scores") {    
        ist.flavours <- c("isd", "isa")
    }    
        
    plyr::ldply(
        stats::setNames(ist.flavours, ist.flavours), 
        function(flavour) {
            mod <- methods::slot(object, flavour)
            
            if (length(mod@genes.mod) > 0) {
                predict(
                    mod, newdata = newdata, 
                    type = type, ncomp = ncomp, ...)
            } else {
                message("Flavour ", flavour, " is not fitted. Skipping...")
            }
        }, 
        .id = "flavour"
    )
}

#' @description \code{predict.ist.pathways()} applies all the fitted 
#' oc/bin models within the object to 
#' new data (in-silico-treated or not)
#' 
#' @rdname predict-ist
#' 
#' @importFrom methods slot
#' @importFrom plyr ldply
#' @importFrom checkmate assertMatrix assertClass qassert
#' @importFrom stats setNames
#' @importFrom BiocParallel bplapply SerialParam
predict.ist.pathways <- function(
    object, 
    newdata, 
    type = "decision",
    ncomp = c(1, 2),
    ...) {
    checkmate::assertMatrix(
        newdata, mode = "numeric", any.missing = FALSE, 
        row.names = "unique", col.names = "unique")
    
    checkmate::assertClass(object, "ist.pathways")
    
    checkmate::qassert(type, "S1")
    checkmate::assertSubset(type, c("decision", "scores"))
    
    checkmate::qassert(ncomp, "N+")
    
    if (type == "decision") {
        ist.mod <- c("list.oc", "list.bin")
    } else if (type == "scores") {    
        ist.mod <- "list.bin"
    }    
    
    # oc, bin
    plyr::ldply(
        stats::setNames(ist.mod, ist.mod), 
        function(mod) {
            list.mod <- methods::slot(object, mod)
            
            # actual models
            ind.mod <- !sapply(list.mod, is.null)
            
            # predict over models
            list.out <- BiocParallel::bplapply(
                list.mod[ind.mod], 
                predict, 
                newdata = newdata, 
                type = type, 
                ncomp = ncomp, 
                BPPARAM = BiocParallel::SerialParam(), 
                ...
            )
            
            plyr::ldply(list.out, identity, .id = "pathway")
        }, .id = "mod"
    )
}
