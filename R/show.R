#' @title Show methods
#' 
#' @description \code{show.ist.translator()} shows an overview of an
#' \code{ist.translator} object.
#' 
#' @param object an \code{ist.translator} or \code{ist.discriminator}
#' 
#' @return \code{NULL}
#' 
#' @name show-ist
#' @rdname show-ist
show.ist.translator <- function(object) {
    
    cat("* ist.translator object\n")
    cat("\n  + Translating from", object@org.from, "to", object@org.to, "\n")
    cat("\n  + Original fold changes\n")
    
    cat("    - Length:", nrow(object@animal.data.table), "\n")
    cat("    - Summary:\n")
    print(summary(object@animal.data.table))

    len.map <- length(object@mapping.data)
    if (len.map == 0) {
        cat("\n  + Data not mapped, or no genes mapping to orthologs.", 
            "Did you fit()?\n")
    } else {
        cat("\n  + Mapped fold changes\n")
        cat("    - Length:", len.map, "\n")
        cat("    - Summary:\n")
        print(summary(object@mapping.data))
    }
    
    invisible()
}

#' @description \code{show.ist.signatures()} shows an overview of an
#' \code{ist.signatures} object.
#' 
#' @rdname show-ist
show.ist.signatures <- function(object) {
    
    cat("* ist.signatures object\n")
    org.sigs <- table(getMetaSig(object)$sig.org)
    
    cat("\n  + Translating", sum(org.sigs), "signature(s) from", 
        paste(names(org.sigs), org.sigs, sep = ":", collapse = ", "))
    cat("\n  + Translating to", object@org.to, "\n")
    
    fc.summaries <- sapply(object@list.sig, function(x) summary(x$logFC))
    if (length(object@list.sig) < 6) {
        cat("  + Fold change summaries:\n")
        print(fc.summaries)
    } else {
        cat("  + Fold change summary of summaries:\n")
        print(apply(fc.summaries, 1, summary))
    }
    
    len.data <- length(object@list.data)
    if (len.data == 0) {
        cat("\n  + Data not mapped, or no genes mapping to orthologs.", 
            "Did you fit()?\n")
    } else {
        cat("\n  + Number of original genes, mapped genes, and coverage\n")

        size.before <- sapply(object@list.sig, nrow)
        size.after <- sapply(object@list.data, length)
        
        size.stats <- list(
            SizeOriginal = size.before, 
            SizeMapped = size.after, 
            MapingRatio = size.after/size.before
        )
        
        print(sapply(size.stats, summary))
    }
    
    invisible()
}

#' @title Show logical ids
#' 
#' @description \code{show.logical.ids()} is a helper to summarise  
#' logical vectors indicating samples
show.logical.ids <- function(object, ids, slot, prefix) {
    if (missing(ids)) ids <- methods::slot(object, slot)
    cat(prefix, 
        "(", 
        sum(ids), 
        "):", 
        head(rownames(object@X)[ids]), 
        "... [y:", 
        head(object@y[ids]), 
        "...]", "\n")
}
    

#' @description \code{show.ist.discriminator()} shows an overview of an
#' \code{ist.discriminator} object.
#' 
#' @rdname show-ist
show.ist.discriminator <- function(object) {
    
    cat("* ist.discriminator object\n")
    cat("\n  + Translating to", object@org.to, "\n")
    cat("\n  + Original data\n")
    cat("    - Dim:", dim(object@X), "\n")
    cat("    - Differential genes (user-defined):", 
        length(object@de.genes), "\n")
    
    # samples
    show.logical.ids(
        object, slot = "id.ref", prefix = "    - Reference samples")
    show.logical.ids(
        object, slot = "id.ist", prefix = "    - Samples to in silico treat")
    show.logical.ids(
        object, slot = "id.oc", prefix = "    - One-class model samples")
    show.logical.ids(
        object, slot = "id.bin", prefix = "    - Binary model samples")
    
    msg.nofit <- "not fitted\n"
    msg.fit <- "fitted\n"
    
    cat("\n  + ISS: ")
    if (length(object@iss@genes.mod) == 0) {
        cat(msg.nofit)
    } else {
        cat(msg.fit, "    - Features:", length(object@iss@genes.top), "\n")
        cat("     - Samples:", length(object@iss@svm.best$fitted), "\n")
    }
    
    cat("\n  + ISD: ")
    if (length(object@isd@genes.mod) == 0) {
        cat(msg.nofit)
    } else {
        cat(msg.fit, "    - Features:", length(object@isd@genes.mod), "\n")
        cat("     - Samples:", dim(object@isd@pls.mod$fitted.values)[1], "\n")
        cat("     - Ncomp:", object@isd@Ncomp, "\n")
    }
    
    cat("\n  + ISA: ")
    if (length(object@isa@genes.mod) == 0) {
        cat(msg.nofit)
    } else {
        cat(msg.fit, "    - Features:", length(object@isa@genes.mod), "\n")
        cat("     - Samples:", dim(object@isa@pls.mod$fitted.values)[1], "\n")
        cat("     - Ncomp:", object@isa@Ncomp, "\n")
    }
    
    invisible()
}

#' @description \code{show.ist.pathways()} shows an overview of an
#' \code{ist.pathways} object.
#' 
#' @rdname show-ist
#' 
#' @importFrom stats quantile
show.ist.pathways <- function(object) {
    
    cat("* ist.pathways object\n")
    cat("\n  + Translating to", object@org.to, "\n")
    cat("\n  + Original data\n")
    cat("    - Dim:", dim(object@X), "\n")
    
    show.logical.ids(
        object, slot = "id.oc", prefix = "    - One-class model samples")
    show.logical.ids(
        object, slot = "id.bin", prefix = "    - Binary model samples")
    
    cat("    - OC pathways:", 
        length(unique(object@pathways.table.oc$path.id)), 
        "pathways,", 
        length(unique(object@pathways.table.oc$gene.id)), 
        "unique genes\n")
    cat("    - BIN pathways:", 
        length(unique(object@pathways.table.bin$path.id)), 
        "pathways,", 
        length(unique(object@pathways.table.bin$gene.id)), 
        "unique genes\n")
    
    msg.nofit <- "not fitted\n"
    msg.fit <- "fitted\n"
    
    cat("\n  + OC: ")
    if (length(object@list.oc) == 0) {
        cat(msg.nofit)
    } else {
        is.oc <- vapply(
            object@list.oc, inherits, FUN.VALUE = TRUE, "ist.mod.oc")
        sv.oc <- vapply(object@list.oc[is.oc], function(x) {
            x@svm.best$tot.nSV
        }, FUN.VALUE = 1L)
        
        cat(msg.fit, "   - Models:", sum(is.oc), "/", length(is.oc), "\n")
        cat("    - Number of support vectors, quartiles (0/1/2/3/4):", 
            paste(
                signif(round(stats::quantile(sv.oc), 3)), 
                collapse = "/"), "\n")
        
        # show(object@iss@svm.best)
    }
    
    cat("\n  + BIN: ")
    if (length(object@list.bin) == 0) {
        cat(msg.nofit)
    } else {
        is.bin <- vapply(
            object@list.bin, inherits, FUN.VALUE = TRUE, "ist.mod.bin")
        ncomp.bin <- vapply(object@list.bin[is.bin], function(x) {
            x@Ncomp
        }, FUN.VALUE = 1)
        
        cat(msg.fit, "   - Models:", sum(is.bin), "/", length(is.bin), "\n")
        cat("    - Number of components, quartiles (0/1/2/3/4):", 
            paste(
                signif(round(stats::quantile(ncomp.bin), 3)), 
                collapse = "/"), "\n")
    }
    
    invisible()
}

#' @description \code{show.ist.results()} shows an overview of an
#' \code{show.ist.results} object.
#' 
#' @rdname show-ist
#' 
#' @importFrom stats quantile
show.ist.results <- function(object) {
    
    cat("* ist.results object\n")
    
    cat("\n  + ist.signatures:", 
        length(object@ist.signatures@list.data), "signatures\n")
    cat("\n  + ist.pathways:", 
        "OC", length(object@ist.pathways@list.oc), 
        ", BIN", length(object@ist.pathways@list.bin), "\n")
    
    show.logical.ids(
        object@ist.pathways, ids = object@id.ref, 
        prefix = "    - Reference samples")
    show.logical.ids(
        object@ist.pathways, ids = object@id.ist, 
        prefix = "    - Samples to in silico treat")
    
    cat("\n  + Decisions OC:", nrow(object@tab.decisions.oc), "rows\n")
    print(head(object@tab.decisions.oc, 3))
    
    cat("\n  + Decisions BIN:", nrow(object@tab.decisions.bin), "rows\n")
    print(head(object@tab.decisions.bin, 3))
    
    cat("\n  + Metrics OC:", nrow(object@tab.metrics.oc), "rows\n")
    print(head(object@tab.metrics.oc, 3))
    
    cat("\n  + Metrics BIN:", nrow(object@tab.metrics.bin), "rows\n")
    print(head(object@tab.metrics.bin, 3))
    
    cat("\n  + PLS weights:", nrow(object@tab.pls.bin), "rows\n")
    print(head(object@tab.pls.bin, 3))
    
    
    invisible()
}