#' @title Fit methods
#' 
#' @description \code{fit.ist.translator()} fits the translator by retrieving
#' orthologs and intersecting the identifiers with those in the fold changes.
#' If the translator was initialised with an ortholog mapping, 
#' it will be reused, this can avoid the slow, repeated querying of the same 
#' orthology mapping when dealing with multiple translators.
#' 
#' @param object an \code{ist.translator} or an 
#' \code{ist.discriminator} object
#' @param ... in \code{fit.ist.discriminator()}, 
#' these parameters are passed to 
#' \code{fit.mod.oc()}, \code{fit.mod.bin()} and optionally 
#' to \code{get.de.genes()}. 
#' Ignored otherwise
#' 
#' @return a fitted object with class \code{ist.translator}, 
#' \code{ist.discriminator}, \code{ist.mod.oc} or \code{ist.mod.bin} 
#' 
#' @name fit-ist
#' @rdname fit-ist
#' 
#' @import data.table
#' @importFrom stats setNames
fit.ist.translator <- function(object, ...) {
    org.from <- object@org.from
    org.to <- object@org.to
    
    # Add mapping table
    if (nrow(object@mapping.data.table) > 0) {
        message("Using predefined orthology mapping")
        tab.map <- data.table::data.table(object@mapping.data.table)
    } else {
        tab.map <- orthIDcon(org.from, org.to)
    }
    
    # Compute vector with (nonzero) fold changes
    tab.fc <- object@animal.data.table
    tab.fc$gene <- rownames(tab.fc)
    
    # only high-confidence orthologs (if confidence is available)
    # confidence tends to be binary, but used .95 
    # just in case
    if ("confidence" %in% colnames(tab.map))
        tab.map <- tab.map[confidence > .95]
    
    tab.join <- tab.map[as.data.table(tab.fc), on = "gene", nomatch = 0]
    
    v.fc <- stats::setNames(tab.join$logFC, tab.join$ortholog)
    
    object@mapping.data.table <- tab.map
    object@mapping.data <- v.fc
    
    if (length(v.fc) == 0) {
        warning("None of the gene ids mapped to the org.from organism")
    } else {
        message(
            "Mapped ", length(v.fc), " fold changes overlaying ", 
            nrow(tab.map), " orthology mappings on ", 
            nrow(tab.fc), " fold changes")
    }
    
    object
}

#' @description \code{fit.ist.signatures()} fits the translator on multiple
#' signatures by first stratifying by organism, retrieving the
#' orthologs and intersecting the identifiers with those in the fold changes.
#' A custom orthology mapping can be used.
#' 
#' @name fit-ist
#' @rdname fit-ist
#' 
#' @import data.table
#' @importFrom stats setNames
#' @importFrom BiocParallel bpmapply
fit.ist.signatures <- function(object, ...) {
    org.from <- getMetaSig(object)$sig.org
    org.to <- object@org.to
    
    object@tab.meta <- data.table::as.data.table(object@tab.meta)
    
    # Add mapping table
    if (length(object@list.mapping) > 0) {
        message("Using predefined orthology mapping")
        list.mapping <- object@list.mapping
    } else {
        message("Computing orthology mappings")
        list.mapping <- lapply(
            unique(org.from), 
            orthIDcon, 
            org.to = org.to
        )
    }
    
    message("Mapping orthologs")
    object@list.data <- BiocParallel::bpmapply(
        function(sig, list.mapping, org) {
            # data.frame with fold changes
            sig$gene <- rownames(sig)
            # orthology mapping
            tab.map <- data.table::as.data.table(list.mapping[[org]])
            
            # only high-confidence orthologs (if confidence is available)
            if ("confidence" %in% colnames(sig))
                sig <- sig[confidence > .95]
            
            tab.join <- tab.map[
                data.table::as.data.table(sig), on = "gene", nomatch = 0]
            
            v.fc <- stats::setNames(tab.join$logFC, tab.join$ortholog)
            
            v.fc
        }, 
        sig = object@list.sig, 
        org = org.from, 
        MoreArgs = list(list.mapping = list.mapping), 
        SIMPLIFY = FALSE,
        ...
    )
    
    object
}

#' @description \code{fit.ist.discriminator()} fits the 
#' discriminator on the org.to
#' organism
#' 
#' @param flavour flavour of in silico treatment; 
#' can contain iss (one-class) and/or 
#' isd, isa (binary classifier)
#' @param ist.translator an \code{ist.translator} object (for isa), 
#' to obtain the list of genes allowed in the model
#' 
#' @rdname fit-ist
#' 
#' @importFrom checkmate qassert assertSubset assertClass
#' @importFrom utils head
fit.ist.discriminator <- function(object, flavour, ist.translator, ...) {
    
    checkmate::qassert(flavour, "S+")
    checkmate::assertSubset(flavour, c("iss", "isd", "isa", "isc"))
    
    # find de.genes if not provided
    if (length(object@de.genes) == 0) {
        message("de.genes not supplied. Computing using limma...")
        object@de.genes <- get.de.genes(object@X, object@y, ...)
    }
    
    if ("iss" %in% flavour) {
        message("Fitting iss...")
        
        X.ref <- object@X[object@id.oc, , drop = FALSE]
        object@iss <- fit.ist.oc(
            X.ref, 
            de.genes = object@de.genes, 
            ...
        )
    } 
    
    if ("isd" %in% flavour) {
        message("Fitting isd...")
        
        object@isd <- fit.ist.bin(
            X = object@X[object@id.bin, , drop = FALSE], 
            y = object@y[object@id.bin, drop = FALSE], 
            de.genes = object@de.genes, 
            ...)
    }
    
    if ("isa" %in% flavour) {
        message("Fitting isa...")
        
        checkmate::assertClass(
            ist.translator, 
            classes = "ist.translator", 
            null.ok = FALSE)
        
        if (object@org.to != ist.translator@org.to) {
            # iss and isd might already be performed, so the object 
            # will contain them even if organism in isa does not match
            warning(
                "The org.to organisms differ between the translator (", 
                ist.translator@org.to, ") and the discriminator (", 
                object@org.to, "). Not adding isa to object...")
            return(object)
        }
        
        v.fc <- ist.translator@mapping.data
        if (length(v.fc) == 0) {
            warning("No mapping data available, did you fit() the translator?")
            return(invisible())
        }
        
        de.genes <- names(v.fc)
        
        object@isa <- fit.ist.bin(
            X = object@X[object@id.bin, , drop = FALSE], 
            y = object@y[object@id.bin, drop = FALSE], 
            de.genes = de.genes, 
            ...)
    }
    
    object
}

#' @title Intersect pathways to specified genes
#' 
#' @param path.tab data.frame with the pathway ids and genes
#' @param id.gene character, genes to include
#' 
#' @return data.frame with intersected genes
#' 
#' @importFrom dplyr filter
map.pathways.to.ref <- function(path.tab, id.gene, prefix) {
    
    if (is.null(path.tab)) return()
    
    id.inpath <- intersect(id.gene, path.tab$gene.id)
    
    path.tab <- dplyr::filter(path.tab, gene.id %in% id.inpath)
    
    n.data <- length(id.gene)
    n.raw <- nrow(path.tab)
    n.inpath <- length(id.inpath)
    n.indata <- nrow(path.tab)
    
    message(
        prefix, " pathways genes after/before intersecting with data: ", 
        n.indata, "/", n.raw, 
        ". Data coverage: ", n.inpath, "/", n.data)
    
    path.tab
}


fit.helper.path.oc <- function(gene.path, object) {
    gene.nopath <- setdiff(colnames(object@X), gene.path)
    
    #try({ #TODO add on.error option to fit.ist.bin to avoid try
        fit.ist.oc(
            X.ref = object@X[object@id.oc, , drop = FALSE], 
            de.genes = gene.nopath
        )
    #})
}

fit.helper.path.bin <- function(gene.path, object) {
    #try({ #TODO add on.error option to fit.ist.bin to avoid try
        fit.ist.bin(
            X = object@X[object@id.bin, , drop = FALSE], 
            y = object@y[object@id.bin, drop = FALSE],
            de.genes = gene.path
        )
    #})
}

#' @description \code{fit.ist.pathways()} fits the 
#' pathway models on the org.to organism (must have pathway data), 
#' by first intersecting the pathways and then fitting the oc/bin models
#' 
#' @rdname fit-ist
#' 
#' @importFrom plyr dlply
#' @importFrom dplyr filter
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom checkmate qassert assertSubset assertClass assertDisjunct
fit.ist.pathways <- function(object, ...) {
    # browser()
    id.genes.ref <- colnames(object@X)
    
    # map the OC pathways
    df.path.oc <- map.pathways.to.ref(
        object@pathways.table.oc, id.genes.ref, prefix = "OC"
    )
    list.path.oc <- split(df.path.oc$gene.id, df.path.oc$path.id)
    
    # map the BIN pathways
    df.path.bin <- map.pathways.to.ref(
        object@pathways.table.bin, id.genes.ref, prefix = "BIN"
    )
    list.path.bin <- split(df.path.bin$gene.id, df.path.bin$path.id)
    
    # metadata, use dummy data.frame if none was provided
    path.oc <- names(list.path.oc)
    path.bin <- names(list.path.bin)
    
    # add pathway metadata
    path.meta.auto <- data.table::data.table(
        path.id = c(path.oc, path.bin), 
        mod = rep(
            c("oc", "bin"), 
            times = c(length(path.oc), length(path.bin))), 
        n.genes = c(
            sapply(list.path.oc, length), 
            sapply(list.path.bin, length))
    )
    if (nrow(object@pathways.meta) == 0) {
        message("No pathway metadata provided. Adding minimal metadata...")
        object@pathways.meta <- path.meta.auto
    } else {
        message("Using provided pathway metadata...")
        checkmate::assertSubset(path.oc, object@pathways.meta$path.id)
        checkmate::assertSubset(path.bin, object@pathways.meta$path.id)
        
        # reserved colnames
        checkmate::assertDisjunct(
            colnames(object@pathways.meta), c("mod", "n.genes"))
        
        # join tables, keep original order
        path.meta.user <- data.table::data.table(
            object@pathways.meta)[path.id %in% c(path.oc, path.bin), ]
        object@pathways.meta <- path.meta.auto[path.meta.user, on = "path.id"]
    }
    
    message("Fitting oc...")
    object@pathways.table.oc <- df.path.oc
    object@list.oc <- BiocParallel::bplapply(
        list.path.oc, 
        fit.helper.path.oc, 
        object = object, 
        BPPARAM = BiocParallel::SerialParam()
    )
    
    message("Fitting bin...")
    object@pathways.table.bin <- df.path.bin
    object@list.bin <- BiocParallel::bplapply(
        list.path.bin,
        fit.helper.path.bin,
        object = object, 
        BPPARAM = BiocParallel::SerialParam()
    )
    
    object
}

#' @description \code{fit.ist.oc()} fits the one-class 
#' support vector machine
#' 
#' @param X.ref a matrix containing the target organism and group for iss
#' @param de.genes disease genes, which will be EXCLUDED from the 
#' training in iss, USED in isd and IGNORED in isa
#' @param n.genes numeric, maximum number of genes to include in the iss model. 
#' Currently genes with the maximum mean expression are kept
#' @param seed integer, seed before train and final model fit, or \code{NA} 
#' not to set any
#' @param gamma,nu,kernel,tunecontrol parameters for iss, passed to 
#' \code{e1071::tune.svm}
#' 
#' @rdname fit-ist
#' 
#' @importFrom e1071 tune.control tune.svm svm
#' @importFrom stats setNames
#' @importFrom utils tail
#' 
#' @export
fit.ist.oc <- function(
    X.ref, 
    de.genes, 
    n.genes = 500, 
    seed = 1,
    args.oc.svm = default.args.oc.svm(), 
    args.oc.train = default.args.oc.train(), 
    ...) {
    
    # non-disesase genes
    genes.mod <- setdiff(colnames(X.ref), de.genes)
    n.iss <- length(genes.mod)
    n.th <- 5
    
    if (n.iss < n.th) {
        warning(
            "Only ", n.iss, " genes (from ", ncol(X.ref), 
            ") are available to fit the one-class model. ", 
            "At least ", n.th, " are required. Returning NULL...")
        return(invisible())
    }
    
    n.min <- 3
    if (nrow(X.ref) < n.min) {
        warning(
            "Only ", nrow(X.ref), " samples are available ", 
            "to fit the one-class model. ", 
            "At least ", n.min, " are required. Returning NULL...")
        return(invisible())
    }
    
    checkmate::qassert(seed, "x1")
    
    genes.top <- colMeans(X.ref[, genes.mod, drop = FALSE]) 
    genes.top <- names(utils::tail(sort(genes.top), n.genes))
    
    X <- as.data.frame(X.ref[, genes.top, drop = FALSE])
    
    if (!is.na(seed)) set.seed(seed)
    
    arg.train <- c(
        list(
            x = X, 
            type = "one-classification", 
            probability = FALSE, # true breaks down 
            tunecontrol = args.oc.train), 
        args.oc.svm)
    
    svm.tuned <- do.call(e1071::tune.svm, arg.train)
    
    arg.best <- c(
        list(x = X, type = "one-classification", kernel = args.oc.svm$kernel), 
        as.list(svm.tuned$best.parameters)
    )
    svm.best <- do.call(e1071::svm, arg.best)
    
    ist.mod.oc(
        genes.mod = genes.mod, 
        genes.top = genes.top,
        id.mod = rownames(X),
        svm.tuned = svm.tuned, 
        svm.best = svm.best)
}

#' @description \code{fit.ist.bin()} fits the binary classifier, 
#' a PLS regression
#' 
#' @param X matrix with the regressors for isd/isa
#' @param y numeric response for isd/isa
#' @param args.bin.train list with options for 
#' \code{caret::train} in isd/isa 
#' @param args.bin.pls list with options for \code{pls::plsr} in isd/isa
#' 
#' @rdname fit-ist
#' 
#' @importFrom caret train trainControl
#' @importFrom pls plsr
#' 
#' @source https://stackoverflow.com/questions/35907477/caret-package-
#' stratified-cross-validation-in-train-function
#' 
#' @export
fit.ist.bin <- function(
    X, 
    y,
    de.genes, 
    seed = 1,
    args.bin.train = default.args.bin.train(),
    args.bin.pls = default.args.bin.pls(), 
    ...) {
    
    genes.mod <- intersect(colnames(X), de.genes)
    
    n.th <- 5
    if (length(genes.mod) <= n.th) {
        warning(
            "5 or less genes in de.genes map to X's colnames (", 
            paste(genes.mod, collapse = ","), 
            "Returning NULL...")
        return(invisible())
    }
    
    n.min <- 3
    if (nrow(X) < n.min) {
        warning(
            "Only ", nrow(X), " samples are available ", 
            "to fit the one-class model. ", 
            "At least ", n.min, " are required. Returning NULL...")
        return(invisible())
    }
    
    checkmate::qassert(args.bin.train$tuneLength, "X1")
    checkmate::qassert(seed, "x1")
    
    # non-disesase genes
    Xde <- X[, genes.mod, drop = FALSE]

    if (!is.na(seed)) set.seed(seed)
    # cross validation to optimise ncomp
    arg.train <- c(
        list(x = Xde, y = y, method = "pls"), 
        args.bin.train)
    # supress the following warning only:
    # You are trying to do regression and your outcome only has two 
    # possible values Are you trying to do classification? 
    # If so, use a 2 level factor as your outcome column.
    withCallingHandlers(
        expr = {pls.train <- do.call(caret::train, arg.train)}, 
        warning = function(w) {
            if (grepl("trying to do classification", w$message))
                invokeRestart("muffleWarning")
        }
    )
    
    # final model
    # 
    # fit all the components specified in the input
    # but predict() will use the optimal number
    arg.mod <- c(
        list(y ~ Xde, ncomp = args.bin.train$tuneLength), 
        args.bin.pls)
    pls.mod <- do.call(pls::plsr, arg.mod)
    
    ist.mod.bin(
        genes.mod = genes.mod, 
        id.mod = rownames(Xde), 
        pls.mod = pls.mod,
        Ncomp = pls.train$finalModel$ncomp, 
        train = pls.train)
    
}

#' @title Compute signatures in long format (after ortholog mapping) 
#' 
#' @param ist.signatures ist.signatures object
#' 
#' @description \code{compute.tab.signatures()} returns the signatures, after 
#' mapping orthologs, in a single data.table
#' 
#' @importFrom data.table data.table
#' @importFrom plyr ldply
#' @importFrom checkmate assertClass
compute.tab.signatures <- function(ist.signatures) {
    checkmate::assertClass(ist.signatures, "ist.signatures")
    
    df.sig <- plyr::ldply(
        ist.signatures@list.data, 
        function(x) data.frame(
            ortholog = names(x), 
            logFC = x, 
            stringsAsFactors = FALSE), 
        .id = "sig.id")
    
    data.table::data.table(df.sig)
}

#' @title Compute the tables with deltas for heatmaps
#' 
#' @param ist.results an ist.results object (needs to have 
#' decisions, pls and signture tables)
#' 
#' @description \code{compute.tab.delta()} returns the table that
#' generats the recovery heatmaps
#' 
#' @importFrom stats median
#' @importFrom checkmate assertClass
compute.tab.delta <- function(ist.results, ...) {
    checkmate::assertClass(ist.results, "ist.results")
    
    sample.names <- rownames(ist.results@ist.pathways@X)
    group.names <- getGroupLevels(ist.results)
    
    # extract necessary table
    df.dec <- getTabDecisions(ist.results)
    df.pls <- getTabPls(ist.results)
    df.sig <- getTabSignatures(ist.results)
    
    # sample identifiers for reference samples and IST'd samples
    nm.ref <- sample.names[ist.results@id.ref]
    nm.ist <- sample.names[ist.results@id.ist]
    
    # take original signature decisions (reference and ISD'd samples)
    df.ref.ist <- df.dec[
        sig.id %in% group.names & sample %in% c(nm.ref, nm.ist)]
    df.ref.ist[
        , type := ifelse(sample %in% nm.ref, "median.ref", "median.ist")]
    
    # compute medians per pathway and group (ref/ist'd)
    df.medians <- df.ref.ist[
        , keyby = .(type, sig.id, pathway), 
        .(decision.median = stats::median(decision.value))]
    
    # cast medians, compute pathway-wise values for 100% recovery
    df.delta100 <- data.table::dcast(
        df.medians, pathway ~ type, value.var = "decision.median")
    df.delta100[, delta.100percent := median.ref - median.ist]
    
    # data.frame for plotting
    lab.na <- "NotSignificant"
    levels.sig <- c(levels(df.sig$sig.id), lab.na)
    
    # must add allow.cartesian because every gene will appear several times 
    # in the signatures, so the joined table can multiply the number of rows 
    # in df.delta100[df.pls], and this is intended. See for instance
    # https://stackoverflow.com/questions/23809517
    df.delta <- df.sig[
        df.delta100[df.pls, on = "pathway", allow.cartesian = TRUE], 
        on = "ortholog", allow.cartesian = TRUE]
    df.delta[, delta := Coef*logFC/Scale]
    # labels for percentages
    df.delta[, delta.percent := delta/delta.100percent*100]
    df.delta[, delta.label := percent.to.labels(delta.percent)]
    # labels for signatures
    df.delta[
        , sig.label := ifelse(is.na(sig.id), lab.na, as.character(sig.id))]
    df.delta[, sig.label := factor(sig.label, levels = levels.sig)] 
    
    # group by pathways, sum contributions
    df.deltapathway <- df.delta[
        !is.na(sig.id), 
        keyby = .(sig.id, pathway), 
        .(total.delta.percent = sum(delta.percent))]
    df.deltapathway[
        , total.delta.label := percent.to.labels(total.delta.percent)]
    
    list(
        tab.delta.median = df.medians, 
        tab.delta.gene = df.delta, 
        tab.delta.pathway = df.deltapathway)
}

#' @description \code{fit.ist.results()} applies a battery of signatures
#' in an \code{ist.signatures} object to the pathway models in 
#' \code{ist.pathways}. Given the complexity of the results, an object of 
#' class \code{ist.results} is returned.
#' 
#' @rdname fit-ist
#' 
#' @import data.table
#' @importFrom checkmate assertClass assertCharacter assertTRUE
#' @importFrom BiocParallel bplapply SerialParam
fit.ist.results <- function(object, ...) {
    # browser()
    ist.path <- object@ist.pathways
    ist.sig <- object@ist.signatures
    
    checkmate::assertClass(ist.sig, "ist.signatures")
    checkmate::assertClass(ist.path, "ist.pathways")
    
    # check there are samples in both ref and ist groups
    checkmate::assertTRUE(sum(object@id.ref) > 0)
    checkmate::assertTRUE(sum(object@id.ist) > 0)
    
    # define sample grops for plos
    if (length(object@group) == 0) {
        message("No groups provided. Grouping by y")
        # must sort them first, otherwise negatives are sorted funny
        # -1, -2, -3, 0, 1, 2, 3.. (should be -3, -2, -1, 0, 1, 2, 3)
        y <- object@ist.pathways@y
        group.lv <- as.character(sort(unique(y)))
        object@group <- ordered(
            paste0("y=", y), 
            levels = paste0("y=", group.lv))
    }
    object@group <- droplevels(object@group)
    
    checkmate::assertCharacter(
        object@vec.gene2label, null.ok = TRUE, unique = TRUE)
    checkmate::assertCharacter(
        names(object@vec.gene2label), null.ok = TRUE, unique = TRUE)
    
    # sigatures
    object@tab.signatures <- compute.tab.signatures(ist.sig)
    
    message("Adding model performance and weights")
    # model tables (weights, perf)
    object@tab.metrics.oc <- data.table::rbindlist(
        lapply(ist.path@list.oc, get.perf.oc), idcol = "pathway")
    object@tab.metrics.bin <- data.table::rbindlist(
        lapply(ist.path@list.bin, get.perf.bin), idcol = "pathway")

    # pls dataframe
    object@tab.pls.bin <- data.table::rbindlist(
        lapply(ist.path@list.bin, get.plstab.bin), idcol = "pathway")
    
    # predict on original (separated by group) and IST'd data
    # convert to df to split, and then back to matrix to predict
    X.orig <- ist.path@X
    list.orig <- split(as.data.frame(X.orig), getGroups(object), drop = TRUE)
    list.X.ist <- c(
        lapply(list.orig, as.matrix), 
        predict(ist.sig, newdata = X.orig[object@id.ist, , drop = FALSE])
    )
    
    message("Computing decision values")
    tab.decisions <- data.table::rbindlist(
        BiocParallel::bplapply(
            list.X.ist, 
            predict, 
            object = ist.path, 
            type = "decision", 
            BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
        ), 
        idcol = "sig.id"
    )
    # asdf
    levels.sig <- c(getGroupLevels(object), getMetaSig(ist.sig)$sig.id)
    tab.decisions[, sig.label := factor(sig.id, levels = levels.sig)]
    
    object@tab.decisions.oc <- tab.decisions[mod == "list.oc"]
    object@tab.decisions.bin <- tab.decisions[mod == "list.bin"]
    
    message("Computing delta tables")
    list.delta <- compute.tab.delta(object)
    object@tab.delta.median <- list.delta$tab.delta.median
    object@tab.delta.gene <- list.delta$tab.delta.gene
    object@tab.delta.pathway <- list.delta$tab.delta.pathway
    
    # add all missing genes to mapping
    genes.background <- sort(unique(list.delta$tab.delta.gene$ortholog))
    genes.nomatch <- setdiff(genes.background, names(object@vec.gene2label))
    labels.nomatch <- paste0("- (", genes.nomatch, ")")
    v.nomatch <- setNames(labels.nomatch, genes.nomatch)
    message("Adding ", length(v.nomatch), " missing genes to vec.gene2label")
    object@vec.gene2label <- c(object@vec.gene2label, v.nomatch)
    
    object
}