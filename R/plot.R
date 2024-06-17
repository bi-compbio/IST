#' @title Plot methods
#' 
#' @description \code{plot.ist.all()} runs on 
#' fitted \code{ist.translator} and \code{ist.discriminator} objects
#' 
#' @param x,y for \code{plot.ist.all()}, respectively fitted 
#' \code{ist.translator} and \code{ist.discriminator} objects
#' @param type character, type of plot: \code{"decision"} (default), 
#' \code{"scores"}
#' @param ... ignored
#' 
#' @return a \code{ggplot} object (plotted as well)
#' 
#' @name plot-ist
#' @rdname plot-ist
#' 
#' @importFrom plyr ldply
#' @importFrom dplyr mutate filter
#' @importFrom checkmate assertClass qassert assertSubset
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
plot.ist.all <- function(x, y, type = "decision", ...) {
    
    checkmate::assertClass(x, "ist.translator")
    checkmate::assertClass(y, "ist.discriminator")
    
    checkmate::qassert(type, "S1")
    checkmate::assertSubset(type, c("decision", "scores"))
    
    # untreated and treated data
    X.untrt <- y@X
    X.trt <- predict(x, newdata = X.untrt[y@id.ist, , drop = FALSE])
    
    # reference samples
    nm.ref <- rownames(X.untrt)[y@id.ref]
    
    # scores or decision?
    if (type == "scores") {
        df.scores <- plyr::ldply(
            list(IST = X.trt, Original = X.untrt),
            predict, 
            object = y, 
            type = "scores",
            .id = "treatment") %>%
            dplyr::mutate(class = ifelse(sample %in% nm.ref, "ref", "other"))
        
        df.ref <- dplyr::filter(
            df.scores, treatment == "None" & class == "ref")
        df.other <- dplyr::filter(df.scores, class == "other")
        
        ggplot2::ggplot() + 
            ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
            ggplot2::geom_vline(xintercept = 0, lty = 2, colour = "gray50") +
            ggplot2::stat_density_2d(
                ggplot2::aes(x = Comp1, y = Comp2, fill = stat(level)), 
                geom = "polygon", data = df.ref) +
            ggplot2::scale_fill_gradientn(
                colours = rev(RColorBrewer::brewer.pal(7, "YlGnBu")), 
                name = "Density") +
            ggplot2::geom_point(
                ggplot2::aes(x = Comp1, y = Comp2, colour = treatment), 
                size = 1.5, pch = 20, data = df.other) +
            ggplot2::scale_colour_brewer(
                palette = "Set1", name = "Treatment") +
            ggplot2::facet_grid(~flavour) + 
            ggplot2::theme_bw() +
            ggplot2::theme(aspect.ratio = 1)
    
    } else if (type == "decision") {
        df.decision <- plyr::ldply(
            list(IST = X.trt, Original = X.untrt),
            predict, 
            object = y, 
            type = "decision",
            .id = "treatment") %>%
            dplyr::mutate(class = ifelse(sample %in% nm.ref, "Ref.", "Other"))
        
        ggplot2::ggplot(df.decision, aes(x = class, y = decision.value)) +
            ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
            ggplot2::geom_boxplot() +
            ggplot2::facet_grid(
                flavour~treatment, scales = "free_x", space = "free_x") + 
            ggplot2::xlab("Class") + 
            ggplot2::ylab("Decision value") + 
            ggplot2::theme_bw() +
            ggplot2::theme(aspect.ratio = 1)
    }
}

#' @description \code{plot.ist.boxplots()} displays the
#' pathway-wise boxplots
#' 
#' @param x,y for \code{plot.ist.boxplots()}, respectively fitted 
#' \code{ist.results} and a character with the pathways to plot
#' @param type character, type of plot (default "ggplot" or "pheatmap")
#' @param vars.meta.sig,vars.meta.path for pheatmap only, 
#' NULL not to display any signature or pathway metadata, or metadata 
#' column names otherwise. Pathways are automatically given the \code{n.genes}
#' metadata, with the number of mapped genes.
#' @param text.size numeric, size of text for heatmaps
#' @param args.pheatmap list with further arguments 
#' to pass to \code{pheatmap()}
#' @param mapping call to \code{aes(...)} to overlay aesthetics in the 
#' boxplots (e.g. \code{aes(fill = org.name)})
#' @param facet_rows,facet_cols for \code{plot.ist.boxplots()}, 
#' facetting variables
#' @param main,main.width plot title (character, or NULL) and maximum width 
#' of title wrapping (integer)
#' @param ... ignored at the moment
#' 
#' @return \code{plot.ist.boxplots()} return a list with the plot
#' object and the plotted data
#' 
#' @name plot-ist
#' @rdname plot-ist
#' 
#' @import ggplot2
#' @importFrom checkmate assertClass qassert assertSubset assertCharacter
#' @importFrom methods show
#' @export
plot.ist.boxplots <- function(
    x, y = head(getPathways(x), 4), 
    sig.ids = c(getSignatures(x), getGroupLevels(x)), mapping = NULL, 
    facet_rows = NULL, facet_cols = ggplot2::vars(pathway), 
    main = NULL, main.width = 50, 
    ...) {
    # browser()
    checkmate::assertClass(x, "ist.results")
    checkmate::qassert(y, "S+")
    # checkmate::assertSubset(y, getPathways(x, mod = "bin"))
    
    checkmate::assertCharacter(
        main, len = 1, null.ok = TRUE, any.missing = FALSE)
    checkmate::qassert(main.width, "N1")
    
    # wrap title
    wrap.main <- ggplot2::label_wrap_gen(width = main.width)
    if (is.null(main)) str.main <- NULL
    else str.main <- head(wrap.main(main), 1)
    
    plot.data <- getTabDecisions(x)[sig.id %in% sig.ids & pathway %in% y]
    if (!is.null(mapping))
        plot.data <- getMetaSig(x)[plot.data, on = "sig.id"]
    
    plot.obj <- ggplot2::ggplot(
        plot.data, 
        ggplot2::aes(x = sig.label, y = decision.value)) +
        ggplot2::geom_hline(yintercept = 0, lty = 1, colour = "gray50") +
        ggplot2::geom_hline(yintercept = -1, lty = 2, colour = "gray70") +
        ggplot2::geom_hline(yintercept = 1, lty = 2, colour = "gray70") +
        ggplot2::geom_boxplot(mapping = mapping, lwd = .2, outlier.size = .5) +
        # scale_fill_manual(values = col.anmod) +
        ggplot2::facet_grid(
            rows = facet_rows, cols = facet_cols, 
            labeller = ggplot2::label_wrap_gen(), 
            scales = "free", space = "free") + 
        ggplot2::xlab("") +
        ggplot2::ylab("Decision value") +
        ggplot2::ggtitle(str.main) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.text.y = ggplot2::element_text(angle = 0))
    
    methods::show(plot.obj)
    
    invisible(list(plot.obj = plot.obj, plot.data = plot.data))
}


#' @description \code{plot.ist.genemaps()} displays the
#' gene heatmaps of a selected pathway
#' 
#' @param x,y for \code{plot.ist.genemaps()}, respectively fitted 
#' \code{ist.results} and a character with the pathway to plot
#' @param max.genes in \code{plot.ist.genemaps()}, numeric, 
#' maximum number of genes to display in genemap
#' (defaults to all). Genes will be prioritised using their sum of squares,
#' useful for large heatmaps.
#' 
#' @return \code{plot.ist.genemaps()} return a list with the plot
#' object and the plotted data
#' 
#' @name plot-ist
#' @rdname plot-ist
#' 
#' @import ggplot2
#' @importFrom checkmate assertClass qassert assertSubset assertCharacter
#' @importFrom pheatmap pheatmap
#' @importFrom methods show
#' @export
plot.ist.genemaps <- function(
    x, y, type = "ggplot", sig.ids = getSignatures(x), max.genes = Inf,
    text.size = 3, vars.meta.sig = NULL, main = y, main.width = 50, 
    args.pheatmap = default.args.pheatmap("genemap"), ...) {
    
    checkmate::assertClass(x, "ist.results")
    checkmate::qassert(y, "S1")
    checkmate::assertSubset(y, getPathways(x, mod = "bin"))
    checkmate::qassert(max.genes, "N1[1,]")
    
    args <- list(...)
    
    checkmate::assertSubset(type, c("ggplot", "pheatmap"))
    
    checkmate::assertCharacter(
        main, len = 1, null.ok = TRUE, any.missing = FALSE)
    checkmate::qassert(main.width, "N1")
    
    # wrap title
    wrap.main <- ggplot2::label_wrap_gen(width = main.width)
    if (is.null(main)) str.main <- NULL
    else str.main <- head(wrap.main(main), 1)
    
    if (type == "ggplot") {
        colours <- colours.lmh()
        plot.data <- get.tab.genemap(
            x, id.path = y, sig.ids = sig.ids, 
            long.only = TRUE, max.genes = max.genes)
        if (is.null(plot.data)) return(invisible())
        
        plot.obj <- ggplot2::ggplot(
            plot.data, 
            ggplot2::aes(
                x = ortholog, y = sig.label, 
                fill = delta.percent, 
                label = delta.label)) +
            ggplot2::geom_tile(width = .8, height = .8) +
            ggplot2::geom_text(size = text.size) +
            ggplot2::xlab("") +
            ggplot2::ylab("") +
            ggplot2::ggtitle(str.main) +
            ggplot2::scale_fill_gradient2(
                low = colours["low"], mid = colours["mid"], 
                high = colours["high"], name = "Recovery") +
            ggplot2::theme_bw() +
            theme.heatmap()
        
        methods::show(plot.obj)
        
    } else if (type == "pheatmap") {
        # browser()
        plot.data <- get.tab.genemap(
            x, id.path = y, sig.ids = sig.ids, long.only = FALSE, 
            vars.meta.sig = vars.meta.sig, max.genes = max.genes)
        if (is.null(plot.data)) return(invisible())
        if (is.null(str.main)) str.main <- NA
        
        allargs.pheatmap <- c(
            list(
                mat = plot.data$data.wide, 
                color = heatmaps.pal(), 
                breaks = sym.breaks(plot.data$data.wide), 
                main = str.main, 
                annotation_row = plot.data$row, 
                annotation_col = plot.data$col[, "Weight", drop = FALSE]
            ), 
            args.pheatmap
        )
        
        plot.obj <- do.call(pheatmap::pheatmap, args = allargs.pheatmap)
    }
    
    invisible(list(plot.obj = plot.obj, plot.data = plot.data))
}

#' @description \code{plot.ist.pathwaymaps()} displays the
#' pathway heatmaps of a selection of pathways
#' 
#' @param x,y for \code{plot.ist.pathwaymaps()}, respectively fitted 
#' \code{ist.results} and a character with multiple pathways to plot
#' 
#' @return \code{plot.ist.pathwaymaps()} return a list with the plot
#' object and the plotted data
#' 
#' @name plot-ist
#' @rdname plot-ist
#' 
#' @import ggplot2
#' @importFrom checkmate assertClass qassert assertSubset assertCharacter
#' @importFrom pheatmap pheatmap
#' @importFrom methods show
#' @export
plot.ist.pathwaymaps <- function(x, y = getPathways(x), 
    type = "ggplot", sig.ids = getSignatures(x), 
    text.size = 3, vars.meta.sig = NULL, vars.meta.path = NULL, 
    args.pheatmap = default.args.pheatmap("pathwaymap"), 
    main = NULL, main.width = 50, ...) {
    
    checkmate::assertClass(x, "ist.results")
    checkmate::qassert(y, "S+")
    checkmate::assertSubset(y, getPathways(x, mod = "bin"))
    
    checkmate::assertSubset(type, c("ggplot", "pheatmap"))
    
    checkmate::assertCharacter(
        main, len = 1, null.ok = TRUE, any.missing = FALSE)
    checkmate::qassert(main.width, "N1")
    
    # wrap title
    wrap.main <- ggplot2::label_wrap_gen(width = main.width)
    if (is.null(main)) str.main <- NULL
    else str.main <- head(wrap.main(main), 1)
    
    if (type == "ggplot") {
        colours <- colours.lmh()
        plot.data <- get.tab.pathwaymap(
            x, id.path = y, sig.ids = sig.ids, long.only = TRUE)
        if (is.null(plot.data)) return(invisible())
        
        plot.obj <- ggplot2::ggplot(
            plot.data, 
            ggplot2::aes(
                x = pathway, y = sig.id, 
                fill = total.delta.percent, 
                label = total.delta.label)) +
            ggplot2::geom_tile() +
            ggplot2::geom_text(size = text.size) +
            ggplot2::xlab("") +
            ggplot2::ylab("") +
            ggplot2::ggtitle(str.main) +
            ggplot2::scale_fill_gradient2(
                low = colours["low"], mid = colours["mid"], 
                high = colours["high"], name = "Recovery") +
            ggplot2::theme_bw() +
            theme.heatmap()
        
        methods::show(plot.obj)
        
    } else if (type == "pheatmap") {
        # browser()
        plot.data <- get.tab.pathwaymap(
            x, id.path = y, sig.ids = sig.ids, long.only = FALSE, 
            vars.meta.sig = vars.meta.sig, vars.meta.path = vars.meta.path)
        if (is.null(plot.data)) return(invisible())
        if (is.null(str.main)) str.main <- NA
        
        allargs.pheatmap <- c(
            list(
                mat = plot.data$data.wide, 
                color = heatmaps.pal(), 
                breaks = sym.breaks(plot.data$data.wide), 
                main = str.main, 
                annotation_row = plot.data$row, 
                annotation_col = plot.data$col
            ), 
            args.pheatmap
        )
        
        plot.obj <- do.call(pheatmap::pheatmap, args = allargs.pheatmap)
    }
    
    invisible(list(plot.obj = plot.obj, plot.data = plot.data))
}
