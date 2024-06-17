# use old classes within our classes
methods::setOldClass("tune")
methods::setOldClass("svm")
methods::setOldClass("mvr")
methods::setOldClass("train")
methods::setOldClass("data.table")

#' One-class model class
#'
#' @slot genes.mod character vector, genes eligible to be in the model
#' @slot genes.top character vector, genes in the model
#' @slot svm.tuned \code{tune} object
#' @slot svm.best \code{svm} object
#' 
#' @rdname ist-classes
#' 
#' @exportClass ist.mod.oc
ist.mod.oc <- methods::setClass(
    "ist.mod.oc",
    slots = list(
        genes.mod = "vector",
        genes.top = "vector",
        id.mod = "vector", 
        svm.tuned = "tune", 
        svm.best = "svm"
    ), 
    validity = function(object) {
        checkmate::qassert(object@genes.mod, "S+")
        checkmate::qassert(object@genes.top, "S+")
        checkmate::qassert(object@id.mod, "S+")
        
        TRUE
    }
)

#' Binary classifier model class
#'
#' @slot genes.mod character vector, genes eligible to be in the model
#' @slot pls.mod \code{mvr} object, from package \code{pls}
#' @slot Ncomp numeric, optimal components for the PLS model
#' @slot train \code{train} object, from package \code{caret}
#' 
#' @rdname ist-classes
#' 
#' @exportClass ist.mod.bin
ist.mod.bin <- methods::setClass(
    "ist.mod.bin",
    slots = list(
        genes.mod = "vector",
        id.mod = "vector", 
        pls.mod = "mvr", 
        Ncomp = "numeric", 
        train = "train"
    ), 
    validity = function(object) {
        checkmate::qassert(object@genes.mod, "S+")
        checkmate::qassert(object@id.mod, "S+")
        checkmate::qassert(object@Ncomp, "X?")
        
        TRUE
    }
)


#' Translator class (fold changes and orthology mapping)
#'
#' @slot animal.data.table data.frame with the animal fold changes
#' @slot org.from,org.to character, organism to translate from and to, 
#' should follow the ENSEMBL notation (hsapiens, mmusculus, etc)
#' @slot mapping.data.table data.frame with the orthology mapping 
#' between both organisms
#' @slot mapping.data numeric named vector with the fold changes to 
#' be added
#' 
#' @rdname ist-classes
#' 
#' @importFrom checkmate assertDataFrame assertNames qassert
#' @exportClass ist.translator
ist.translator <- methods::setClass(
    "ist.translator",
    slots = list(
        animal.data.table = "data.frame",
        org.from = "character", 
        org.to = "character", 
        mapping.data.table = "data.frame", 
        mapping.data = "vector"
    ),
    validity = function(object) {
        checkmate::assertDataFrame(
            object@animal.data.table, types = c("numeric"), 
            any.missing = FALSE, row.names = "named", min.rows = 1, 
            col.names = "named", min.cols = 2)
        checkmate::assertNames(
            names(object@animal.data.table), 
            must.include = c("logFC", "adj.P.Val"), 
            disjunct.from = c("confidence", "type"))
        
        checkmate::qassert(object@org.from, "S1")
        checkmate::qassert(object@org.to, "S1")
        
        TRUE
    }
)

#' Signatures class (fold changes, metadata, orthology mapping)
#'
#' @slot list.sig list with data.frames containing the animal fold changes
#' @slot tab.meta data.frame with the metadata - must contain the columns 
#' \code{c("sig.id", "sig.name", "sig.org")}.
#' @slot org.to character, organism to translate to, 
#' should follow the ENSEMBL notation (hsapiens, mmusculus, etc)
#' @slot list.mapping list with data.frames, by organism 
#' (to translate from) name
#' @slot list.data list with numeric named vectors with the fold changes to 
#' be added
#' 
#' @rdname ist-classes
#' 
#' @importFrom checkmate assertDataFrame assertNames assertCharacter 
#' assertSetEqual qassert
#' @exportClass ist.signatures
ist.signatures <- methods::setClass(
    "ist.signatures",
    slots = list(
        list.sig = "list",
        tab.meta = "data.frame", 
        org.to = "character", 
        list.mapping = "list", 
        list.data = "list"
    ),
    validity = function(object) {
        # check signatures
        checkmate::assertCharacter(names(object@list.sig), unique = TRUE)
        
        lapply(
            object@list.sig, checkmate::assertDataFrame, types = "numeric", 
            any.missing = FALSE, min.rows = 1, min.cols = 2, 
            row.names = "named", col.names = "named"
        )
        
        lapply(
            object@list.sig, 
            function(x, ...) checkmate::assertNames(names(x), ...), 
            must.include = c("logFC", "adj.P.Val"), 
            disjunct.from = c("confidence", "type")
        )
        
        # check metadata
        col.meta <- c("sig.id", "sig.name", "sig.org")
        checkmate::assertNames(
            names(object@tab.meta), must.include = col.meta
        )
        
        checkmate::assertDataFrame(
            object@tab.meta[col.meta], nrows = length(object@list.sig), 
            any.missing = FALSE, types = "character"
        )
        
        checkmate::assertSetEqual(
            names(object@list.sig), object@tab.meta$sig.id, ordered = TRUE
        )
        
        # check organism
        checkmate::qassert(object@org.to, "S1")
        
        TRUE
    }
)

#' Discriminator class (one-class/binary model and pathways)
#'
#' @slot X numeric matrix, independent variables
#' @slot y numeric vector, dependent variable
#' @slot org.to character, organism with ENSEMBL notation 
#' (hsapiens, mmusculus, etc)
#' @slot de.genes character, disease-associated genes 
#' (must be a subset of X's colnames)
#' @slot id.oc,id.bin,id.ref,id.ist logical vectors, with length equal 
#' to \code{nrow(X)}, indicating which samples should be used to fit 
#' the one class, binary models, and be treated as reference and target
#' (in silico treatment should ideally move the reference to the target group)
#' @slot pathways.data.table data.frame with pathway annotations
#' @slot iss,isd,isa fitted models
#' 
#' @rdname ist-classes
#' 
#' @importFrom checkmate qassert assertMatrix assertSubset assertDataFrame
#' assertLogical
#' @exportClass ist.discriminator
ist.discriminator <- methods::setClass(
    "ist.discriminator",
    slots = list(
        X = "matrix",
        y = "vector",
        org.to = "character", 
        de.genes = "character",
        id.oc = "logical",
        id.bin = "logical", 
        id.ref = "logical", 
        id.ist = "logical", 
        pathways.data.table = "data.frame", 
        iss = "ist.mod.oc",
        isd = "ist.mod.bin",
        isa = "ist.mod.bin"
    ), 
    validity = function(object) {
        checkmate::assertMatrix(
            object@X, mode = "numeric", any.missing = FALSE, 
            row.names = "unique", col.names = "unique")
        checkmate::qassert(object@y, "N+")
        if (nrow(object@X) != length(object@y)) 
            stop("X must have as many rows as elements in y")
        
        checkmate::qassert(object@org.to, "S1")
        
        checkmate::qassert(object@de.genes, "S*")
        # checkmate::qassert(object@de.genes, min.len = 5)
        checkmate::assertSubset(object@de.genes, colnames(object@X))
        
        qrule.bool <- paste0("B", nrow(object@X))
        checkmate::qassert(object@id.oc, qrule.bool)
        checkmate::qassert(object@id.bin, qrule.bool)
        checkmate::qassert(object@id.ref, qrule.bool)
        checkmate::qassert(object@id.ist, qrule.bool)
        
        checkmate::assertDataFrame(
            object@pathways.data.table, types = "character", 
            any.missing = FALSE, min.rows = 0, col.names = "named", 
            max.cols = 2)
        
        TRUE
    }
)

#' Discriminator class (one-class/binary model and pathways)
#'
#' @slot X numeric matrix, independent variables
#' @slot y numeric vector, dependent variable
#' @slot org.to character, organism with ENSEMBL notation 
#' (hsapiens, mmusculus, etc)
#' @slot pathways.table.oc,pathways.table.bin data.frame with 
#' pathway annotations for pathways expected to be unchanged/changed
#' @slot list.oc,list.bin list with fitted models for 
#' unchanged/changed pathways
#' 
#' @rdname ist-classes
#' 
#' @importFrom checkmate assertDataFrame assertMatrix 
#' assertArray qassert assertList assertSubset assertCharacter
#' @exportClass ist.pathways
ist.pathways <- methods::setClass(
    "ist.pathways",
    slots = list(
        X = "matrix",
        y = "vector",
        org.to = "character", 
        id.oc = "logical",
        id.bin = "logical", 
        pathways.table.oc = "data.frame", 
        pathways.table.bin = "data.frame", 
        pathways.meta = "data.frame", 
        list.oc = "list",
        list.bin = "list"
    ), 
    prototype = list(
        X = matrix(nrow = 0, ncol = 0),
        y = numeric(),
        org.to = character(), 
        id.oc = logical(),
        id.bin = logical(), 
        pathways.table.oc = data.frame(
            path.id = character(), gene.id = character(), 
            stringsAsFactors = FALSE), 
        pathways.table.bin = data.frame(
            path.id = character(), gene.id = character(), 
            stringsAsFactors = FALSE), 
        pathways.meta = data.frame(
            path.id = character(), stringsAsFactors = FALSE), 
        list.oc = list(),
        list.bin = list()
    ),
    validity = function(object) {
        checkmate::assertMatrix(
            object@X, mode = "numeric", any.missing = FALSE, 
            row.names = "unique", col.names = "unique")
        checkmate::assertNumeric(
            object@y, any.missing = FALSE, finite = TRUE, len = nrow(object@X))
        
        checkmate::qassert(object@org.to, "S1")
        
        qrule.bool <- paste0("B", nrow(object@X))
        checkmate::qassert(object@id.oc, qrule.bool)
        checkmate::qassert(object@id.bin, qrule.bool)

        checkmate::assertDataFrame(
            object@pathways.table.oc, types = "character", 
            any.missing = FALSE, min.rows = 0, col.names = "named", 
            max.cols = 2)
        checkmate::assertDataFrame(
            object@pathways.table.bin, types = "character", 
            any.missing = FALSE, min.rows = 0, col.names = "named", 
            max.cols = 2)
        
        # check that at least we have the path.id column in meta
        # unique values, and containing all those in the pathways
        checkmate::assertDataFrame(
            object@pathways.meta, col.names = "named", min.col = 1)
        checkmate::assertSubset("path.id", colnames(object@pathways.meta))
        checkmate::assertCharacter(
            object@pathways.meta$path.id, unique = TRUE, any.missing = FALSE)
        
        
        checkmate::assertList(
            object@list.oc, types = c("ist.mod.oc", "try-error"))
        checkmate::assertList(
            object@list.bin, types = c("ist.mod.bin", "try-error"))
        
        TRUE
    }
)

#' Results class (one-class/binary pathways)
#'
#' @slot ist.signatures,ist.pathways objects to generate the results from
#' @slot vec.gene2label nicer labels for plotting genes
#' @slot group factor, groups to stratify samples in plotting
#' @slot tab.decisions.oc,tab.metrics.oc data.table objects for OC models
#' @slot tab.decisions.bin,tab.metrics.bin,tab.pls.bin data.table objects 
#' for BIN models
#' 
#' @rdname ist-classes
#' 
#' @importFrom data.table data.table
#' @importFrom checkmate assertClass assertFactor assertTRUE qassert
#' @exportClass ist.results
ist.results <- methods::setClass(
    "ist.results",
    slots = list(
        ist.signatures = "ist.signatures",
        ist.pathways = "ist.pathways",
        id.ref = "logical", 
        id.ist = "logical", 
        group = "factor", 
        vec.gene2label = "character", 
        tab.signatures = "data.table", 
        tab.decisions.oc = "data.table", 
        tab.decisions.bin = "data.table", 
        tab.metrics.oc = "data.table", 
        tab.metrics.bin = "data.table", 
        tab.pls.bin = "data.table", 
        tab.delta.median = "data.table",
        tab.delta.gene = "data.table",
        tab.delta.pathway = "data.table"
    ), 
    prototype = list(
        ist.signatures = ist.signatures(),
        ist.pathways = ist.pathways(),
        id.ref = logical(), 
        id.ist = logical(), 
        group = factor(), 
        vec.gene2label = character(), 
        tab.signatures = data.table::data.table(), 
        tab.decisions.oc = data.table::data.table(), 
        tab.decisions.bin = data.table::data.table(), 
        tab.metrics.oc = data.table::data.table(), 
        tab.metrics.bin = data.table::data.table(), 
        tab.pls.bin = data.table::data.table(), 
        tab.delta.median = data.table::data.table(), 
        tab.delta.gene = data.table::data.table(), 
        tab.delta.pathway = data.table::data.table()
    ), 
    validity = function(object) {
        checkmate::assertClass(
            object@ist.signatures, classes = "ist.signatures")
        checkmate::assertClass(
            object@ist.pathways, classes = "ist.pathways")
        
        checkmate::assertFactor(
            object@group, any.missing = FALSE) #, len = nrow(object@X)
        checkmate::assertTRUE(
            length(object@group) %in% c(0, nrow(object@ist.pathways@X)))
        
        qrule.bool <- paste0("B", nrow(object@ist.pathways@X))
        checkmate::qassert(object@id.ref, qrule.bool)
        checkmate::qassert(object@id.ist, qrule.bool)
        
        TRUE
    }
)
