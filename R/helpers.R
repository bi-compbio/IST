#' @title Helper functions
#' 
#' @description \code{define()} allows making \code{new()} objects with 
#' ellipsis, by matching the class slots to the arguments in \code{...} 
#' and ignoring the rest (unlike \code{new()}, which raises error)
#' 
#' @param Class character, as in \code{new()}
#' 
#' @examples 
#' data(sample.data.ist)
#' ## This gives an error because "a" is not a slot in "ist.discriminator"
#' \dontrun{
#' ist.discr <- with(sample.data.ist, {
#' methods::new(
#' "ist.discriminator", 
#' X = X.hsa, 
#' y = y.hsa, 
#' org.to = "hsapiens", 
#' a = 3
#' )
#' })
#' }
#' ## This works
#' ist.discr <- with(sample.data.ist, {
#' define(
#' "ist.discriminator", 
#' X = X.hsa, 
#' y = y.hsa, 
#' org.to = "hsapiens", 
#' id.oc = rep(TRUE, nrow(X.hsa)), 
#' id.bin = rep(TRUE, nrow(X.hsa)), 
#' id.ref = rep(TRUE, nrow(X.hsa)), 
#' id.ist = rep(TRUE, nrow(X.hsa)), 
#' a = 3
#' )
#' })
#' ist.discr
#' 
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom methods new slotNames
#' @importFrom checkmate qassert 
#' @export
define <- function(Class, ...) {
    checkmate::qassert(Class, "S1")
    args <- list(...)
    
    # match slot names with those of class
    match.slots <- intersect(
        names(args), 
        methods::slotNames(Class)
    )
    
    # call new with those slots only
    do.call(methods::new, c(list(Class = Class), args[match.slots]))
}


#' @description \code{reduce.reg()} applies the desired reduce function
#' to a registry with \code{batchtools}
#' 
#' @param reg.name character, name of the registry
#' @param f function to reduce results, of the form \code{f(aggr, res)}, 
#' that specifies how to combine the partially reduced results \code{aggr}
#' with the next result \code{res}. The default value just binds the 
#' rows of every (decisions) data.frame in the registry.
#' 
#' @return \code{reduce.reg()} returns the reduced results, typically
#' a data.frame
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom batchtools loadRegistry reduceResults waitForJobs
#' @importFrom checkmate qassert check_function
#' @export
reduce.reg <- function(
    reg.name, 
    waitForJobs = TRUE,
    f = function(aggr, res) rbind(aggr, res$decisions), 
    ...) {
    
    checkmate::qassert(reg.name, "S1")
    checkmate::qassert(waitForJobs, "B1")
    checkmate::check_function(f, args = c("aggr", "res"))
    
    reg <- batchtools::loadRegistry(reg.name, ...)
    
    if (waitForJobs) {
        message("Waiting for jobs to finish...")
        batchtools::waitForJobs(reg = reg)
    }

    batchtools::reduceResults(fun = f, init = list(), reg = reg)
}

#' @description \code{find.de.genes()} obtains a list of differential 
#' genes from a matrix of expression values and a response 
#' using \code{limma}
#' 
#' @param X,y matrix (genes in columns) and vector with response variable, 
#' which can be quantitative or two-level factor character
#' @param n.min integer, minimum genes to report (not less than 10)
#' @param return.toptable logical, return the topTable or just the 
#' differential genes? (default: \code{FALSE})
#' @param robust logical, whether to use 
#' robust estimates in \code{limma::eBayes} (if \code{TRUE}, the 
#' \code{statmod} package needs to be installed)
#' @param logFC,adj.P.Val numeric, cutoffs for log2FC and adjusted p-value
#' @param ... \code{define()}: fields to populate the new object, as 
#' in \code{new()}; ignored in \code{get.de.genes}; passed to 
#' \code{loadRegistry} in \code{reduce.reg}
#' 
#' @return 
#' \code{find.de.genes()}: if \code{return.toptable = TRUE}, 
#' a list with two slots: the data.frame \code{topTable} 
#' with the whole DE output, 
#' and a character vector \code{de.genes} with the significant IDs.
#' If \code{return.toptable = FALSE}, just a character vector 
#' of \code{de.genes}.
#' 
#' @examples 
#' data(sample.data.ist)
#' ## don't show the renamed intercept warning
#' suppressWarnings({
#' de.genes <- with(sample.data.ist, {
#' get.de.genes(X.hsa, y.hsa)
#' })
#' de.genes.toptable <- with(sample.data.ist, {
#' get.de.genes(X.hsa, y.hsa, return.toptable = TRUE)
#' }) 
#' })
#' str(de.genes)
#' str(de.genes.toptable)
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom utils head
#' @importFrom stats model.matrix
#' @importFrom checkmate qassert assert_matrix assert_vector assert_true
#' @importFrom limma lmFit makeContrasts contrasts.fit topTable eBayes
#' @export
get.de.genes <- function(
    X, y, n.min = 100, return.toptable = FALSE, 
    robust = FALSE, logFC = .25, adj.P.Val = .05, ...) {
    
    checkmate::assert_matrix(X, "numeric", all.missing = FALSE)
    checkmate::assert_vector(y, all.missing = FALSE)
    checkmate::assert_true(inherits(y, c("integer", "numeric", "factor")))
    
    checkmate::qassert(return.toptable, "B1")
    checkmate::qassert(n.min, "N1[10,)")
    checkmate::qassert(robust, "B1")
    checkmate::qassert(logFC, "N1")
    checkmate::qassert(adj.P.Val, "N1")
    
    if (is.factor(y) & nlevels(y) != 2) 
        stop("The response y is a factor but with ", 
        nlevels(y), " levels instead of 2")
    
    design <- stats::model.matrix(~y, data = data.frame(y = y))
    
    fit <- limma::lmFit(t(X), design)
    
    cnt <- limma::makeContrasts(
        contrasts = colnames(design)[2], 
        levels = design)
    
    fit.cnt <- limma::contrasts.fit(fit, cnt)
    ebayes <- limma::eBayes(fit.cnt, robust = robust)
    
    # whole table
    topTable <- limma::topTable(
        ebayes, number = Inf, adjust.method = "fdr")
    
    # significant hits only
    topDE <- limma::topTable(
        ebayes, number = Inf,  adjust.method = "fdr", 
        p.value = adj.P.Val, lfc = logFC)
    
    # number of genes to report
    n.de <- nrow(topDE)
    message("Found ", n.de, " differential genes.")
    if (n.de < n.min) {
        message(
            "Less than ", n.min, " genes are differential. ", 
            "Using top ", n.min, " instead...")
        
        n.de <- n.min
    }
    
    # take top n.de genes from topTable
    de.genes <- rownames(utils::head(topTable, n.de))
    
    if (return.toptable) {
        # return both for traceability
        return(list(topTable = topTable, de.genes = de.genes))
    }
    return(de.genes)
}

#' @title Generate readable labels
#' 
#' @param percentage numeric vector to translate to plot labels
#' @param label.na character to display when percentage is NA
#' 
#' @examples
#' percent.to.labels(c(107.2223, 24.4, NA))
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom methods new slotNames
#' @importFrom checkmate qassert 
#' @export
percent.to.labels <- function(percent, label.na = "-") {
    ifelse(is.na(percent), label.na, paste0(round(percent), "%"))
}

#' @title Colour low/mid/high for heatmaps
#' 
#' @examples
#' colours.lmh()
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @export
colours.lmh <- function() {
    c(low = "deepskyblue3", mid = "white", high = "darkorange1")
}

#' @title Colour palette for heatmaps
#' 
#' @param n numeric, number of colours to return
#' 
#' @examples
#' heatmaps.pal(13)
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom checkmate qassert 
#' @export
heatmaps.pal <- function(n = 19) {
    checkmate::qassert(n, "N1")
    
    grDevices::colorRampPalette(colours.lmh())(n)
}

#' @title Data breaks for heatmaps (symmetric)
#' 
#' @param mat numeric matrix to compute the range from
#' @param n numeric, number of colours in heatmap
#' 
#' @examples
#' sym.breaks(matrix(1:10, nrow = 2))
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom checkmate qassert assertMatrix
#' @export
sym.breaks <- function(mat, n = 19) {
    checkmate::assertMatrix(mat, mode = "numeric", any.missing = FALSE)
    checkmate::qassert(n, "N1")
    
    maxabs <- max(abs(mat))
    seq(-maxabs, maxabs, length.out = n + 1)
}

#' @title Rank signatures or pathways using the pathwaymaps data
#' 
#' @description Metrics are ranked according to three criteria:
#' \describe{
#'   \item{1}{Mean absolute error to 100%}
#'   \item{2}{Mean absolute error to 100%, but not penalising overshoot}
#'   \item{3}{Mean percentage}
#' }
#' A rank is computed for each one, and then the ranks are averaged into 
#' a final rank.
#' Z-scores are also provided.
#' 
#' @param object the \code{ist.results} object holding the results
#' @param what character, what to rank; either \code{"signatures"} 
#' or \code{"pathways"}
#' @param max.out numeric, maximum number of entities to return 
#' (defaults to 5, can be set to \code{Inf}) 
#' @param return.table logical, whether to return the table with 
#' the metrics used for prioritising 
#' @param id.path,sig.ids ids of pathways and signatures to prioritise 
#' (defaults to all)
#' 
#' @return A character vector if \code{return.table = FALSE}, a 
#' \code{data.table} otherwise
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom stats sd median
#' @importFrom caret MAE
#' @importFrom checkmate qassert assertSubset
#' @import data.table
#' @export
rankPathwaymaps <- function(
    object, what = "signatures", 
    max.out = 5, return.table = FALSE,   
    id.path = getPathways(object), 
    sig.ids = getSignatures(object)) {
    
    vec.name2rank <- c(signatures = "sig.id", pathways = "pathway")
    
    checkmate::qassert(what, "S1")
    checkmate::assertSubset(what, names(vec.name2rank))
    
    checkmate::qassert(max.out, "N1")
    checkmate::qassert(return.table, "B1")
    
    var.rankby <- vec.name2rank[what]
    
    # retrieve long table with overall recovery
    dt.long <- get.tab.pathwaymap(
        object, id.path = id.path, 
        sig.ids = sig.ids, long.only = TRUE)
    
    dt.rank <- dt.long[
        , by = var.rankby, .(
            perc_mae100 = caret::MAE(total.delta.percent, 100), 
            perc_mae100trim = caret::MAE(pmin(total.delta.percent, 100), 100), 
            perc_mean = mean(total.delta.percent), 
            perc_median = stats::median(total.delta.percent), 
            perc_sd = stats::sd(total.delta.percent)
        )
        ]
    
    dt.rank[
        , `:=`( 
            rank_mae100 = rank(perc_mae100), 
            rank_mae100trim = rank(perc_mae100trim), 
            rank_mean = rank(-perc_mean),
            z_mae100 = -scale(perc_mae100), 
            z_mae100trim = -scale(perc_mae100trim), 
            z_mean = scale(perc_mean)
        )
        ]
    
    dt.rank[
        , `:=`(
            rank_aggr = (rank_mae100 + rank_mae100trim + rank_mean)/3, 
            z_aggr = (z_mae100 + z_mae100trim + z_mean)/3
        )
        ]
    
    data.table::setorder(dt.rank, rank_aggr)
    
    dt.rank <- head(dt.rank, max.out)
    
    if (return.table) return(dt.rank)
    # otherwise, just the ids
    as.character(dt.rank[[var.rankby]])
}

#' @title ggplot2 theme for heatmaps
#' 
#' @examples
#' theme.heatmap()
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom ggplot2 theme element_blank element_text
#' @export
theme.heatmap <- function() {
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(), 
        axis.text.x = ggplot2::element_text(
            angle = 45, vjust = 1, hjust = 1))
}

#' @title Save a pheatmap object
#' 
#' @param x pheatmap object, or a list with a pheatmap object, 
#' as returned by the plotting functions
#' @param filename output file (with extension)
#' @param width,height numeric, dev size in inches by default
#' @param dev function, graphical device to use
#' @param dev.args list of arguments to the \code{dev} function
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom checkmate assertClass qassert
#' @importFrom grDevices dev.off png
#' @importFrom grid grid.newpage grid.draw
#' @export
save.pheatmap <- function(
    x, filename, dev = grDevices::png, width = 7, height = 7,
    dev.args = list(res = 150, units = "in")) {
    
    if (inherits(x, "list")) {
        message("List provided in x. Using x <- x$plot.obj")
        x <- x$plot.obj
    }
    
    checkmate::assertClass(x, "pheatmap")
    checkmate::qassert(filename, "S1")
    checkmate::qassert(width, "N1")
    checkmate::qassert(height, "N1")
    checkmate::assertClass(dev, "function")
    checkmate::assertClass(dev.args, "list")
    
    dev.args <- c(
        dev.args, 
        list(filename = filename, width = width, height = height))
    
    do.call(dev, dev.args)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    grDevices::dev.off()
}

#' @title Save all the genemaps as pheatmaps
#' 
#' @param ist.results object to plot
#' @param dir.out character, directory to save plots (will be created
#' if it does not exist)
#' @param id.paths character, pathway ids to plot
#' (defaults to all bin pathways)
#' @param cex.width,cex.height numeric values, scaling factor with respect
#' to default dimensions. The defaults are guesses that attempt to
#' display a reasonable plot.
#' @param ... in \code{save.genemaps()}, graphical parameters to 
#' pass to \code{plot.ist.genemaps()}, like \code{sig.ids}
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom stats setNames
#' @importFrom checkmate assertClass qassert
#' @importFrom grDevices dev.off png
#' @importFrom grid grid.newpage grid.draw
#' @export
save.genemaps <- function(
    ist.results, dir.out, 
    id.paths = getPathways(ist.results), 
    cex.width = 1, cex.height = 1, ...) {
    
    checkmate::qassert(dir.out, "S1")
    checkmate::qassert(id.paths, "S+")
    checkmate::qassert(cex.width, "N1")
    checkmate::qassert(cex.height, "N1")
    
    if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)
    
    # generate plot objects
    list.genemap <- lapply(
        stats::setNames(id.paths, id.paths), 
        plot.ist.genemaps, 
        x = ist.results, 
        type = "pheatmap", 
        ...
    )
    
    # take successful plots
    list.genemap <- list.genemap[!sapply(list.genemap, is.null)]
    if (length(list.genemap) == 0) return(invisible())
    
    # define dimensions based on number of signatures and genes in pathway
    list.size <- lapply(list.genemap, function(x) dim(x$plot.data$data.wide))
    v.pathsize <- sapply(list.size, tail, 1)
    v.sigsize <- sapply(list.size, head, 1)
    
    v.files <- paste0(
        dir.out, "/", gsub(" ", "_", names(list.genemap)), ".png")
    
    width.min <- 7
    height.min <- 3
    size.max <- 50
    
    # default dimensions
    v.width <- pmax(pmin(20*v.pathsize/100*cex.width, size.max), width.min)
    v.height <- pmax(pmin(3 + v.sigsize/5*cex.height, size.max), height.min)
    
    # save the png files
    mapply(
        save.pheatmap, 
        x = list.genemap, 
        filename = v.files, 
        width = v.width,
        height = v.height
    )
    
    invisible(list.genemap)
}

#' @title Save boxplots by pages
#' 
#' @description Divide all the specified pathways in grous of \code{k} 
#' pathways and draw their boxplots
#' 
#' @param k integer, number of plots per page
#' @param width,height numeric, figure size (inches)
#' @param ... in \code{save.boxplots()}, parameters passed to 
#' \code{plot.ist.boxplots()}
#' 
#' @name helpers
#' @rdname helpers
#' 
#' @importFrom checkmate qassert
#' @importFrom ggplot2 ggsave
#' @export
save.boxplots <- function(
    ist.results, dir.out, id.paths = getPathways(ist.results), 
    k = 4, width = 8, height = 4, ...) {
    
    checkmate::qassert(dir.out, "S1")
    checkmate::qassert(id.paths, "S+")
    checkmate::qassert(k, "N1[1,)")
    checkmate::qassert(width, "N1")
    checkmate::qassert(height, "N1")
    
    if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)
    
    # split pathways
    n.pages <- ceiling(length(id.paths)/k)
    suppressWarnings(
        list.paths <- split(
            id.paths, 
            paste0("boxplot-", rep(seq_len(n.pages), each = k))
        ) 
    )
    
    # generate ggplot objects
    list.plt <- lapply(
        list.paths, 
        plot.ist.boxplots, 
        x = ist.results, 
        ...
    )
    
    # save them
    lapply(
        names(list.paths), 
        function(x) ggsave(
            filename = paste0(dir.out, "/", x, ".png"), 
            plot = list.plt[[x]]$plot.obj, 
            width = width, height = height
        )
    )
    
    invisible()
}