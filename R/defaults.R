#' @title Default arguments
#' 
#' @description \code{default.args.bin.train()} generates proper default 
#' arguments for \code{caret::train}, also generating an array of seeds 
#' given the amount of folds, repetitions and parameter values
#' 
#' @param tuneLength,number,repeats numeric arguments for \code{caret::train}
#' 
#' @return \code{default.args.*} return a \code{list} with 
#' default arguments for \code{caret::train} and \code{pls::plsr}
#' 
#' @examples 
#' IST::default.args.bin.train()
#' 
#' @name defaults
#' @rdname defaults
#' 
#' @importFrom caret trainControl 
#' @export
default.args.bin.train <- function(
    tuneLength = 5, 
    number = 5, 
    repeats = 20) {
    
    nseeds <- tuneLength*number*repeats
    matseeds <- matrix(seq_len(nseeds), ncol = number*repeats)
    listseeds <- as.list(as.data.frame(matseeds))
    argseeds <- c(listseeds, list(extra = nseeds + 1))
    
    list(
        preProcess = c("center", "scale"),
        tuneLength = tuneLength,
        trControl = caret::trainControl(
            method = "repeatedcv", 
            number = number, 
            repeats = repeats, 
            seeds = argseeds)
    )
}

#' @description \code{default.args.bin.pls()} generates proper default 
#' arguments for \code{pls::plsr}
#' 
#' @param validation,segments,segment.type,jackknife arguments for 
#' \code{pls::plsr}
#' 
#' @examples 
#' IST::default.args.bin.pls()
#' 
#' @name defaults
#' @rdname defaults
#' 
#' @export
default.args.bin.pls <- function(
    validation = "CV",
    segments = 5,
    segment.type = "random",
    jackknife = FALSE) {
    
    list(
        center = TRUE,
        scale = TRUE,
        validation = validation,
        segments = segments,
        segment.type = segment.type,
        jackknife = jackknife
    )
}

#' @description \code{default.args.oc.train()} generates proper default 
#' arguments for \code{e1071::tune.svm}
#' 
#' @param sampling,nrepeat,cross,error.fun arguments for 
#' \code{e1071::tune.control}
#' 
#' @examples 
#' IST::default.args.oc.train()
#' 
#' @name defaults
#' @rdname defaults
#' 
#' @export
default.args.oc.train <- function(
    sampling = "cross",
    nrepeat = 5,
    cross = 5, 
    error.fun = function(true, pred) mean(!pred)) {
    
    e1071::tune.control(
        sampling = sampling,
        nrepeat = nrepeat,
        cross = cross, 
        error.fun = error.fun)
}

#' @description \code{default.args.oc.svm()} generates proper default 
#' arguments for \code{e1071::tune.svm}
#' 
#' @param gamma,nu,kernel arguments for 
#' \code{e1071::tune.svm} and \code{e1071::svm}
#' 
#' @examples 
#' IST::default.args.oc.svm()
#' 
#' @name defaults
#' @rdname defaults
#' 
#' @export
default.args.oc.svm <- function(
    gamma = 2^(c(-10, -5, -1, 0)),
    nu = c(.1, .3, .5, .7, .9),
    kernel = "radial") {
    
    list(gamma = gamma, nu = nu, kernel = kernel)
}

#' @description \code{default.args.pheatmap()} returns proper default 
#' arguments for plotting pheatmap heatmaps
#' 
#' @param plottype character, type of plot to return defaults for 
#' (genemap or pathwaymap)
#' 
#' @examples 
#' IST::default.args.pheatmap("genemap")
#' 
#' @name defaults
#' @rdname defaults
#' 
#' @export
default.args.pheatmap <- function(plottype) {
    if (plottype == "genemap"){
        return(
            list(cluster_rows = FALSE, cluster_cols = FALSE)
        )
    }
        
    
    if (plottype == "pathwaymap") {
        return(
            list(
                cluster_rows = FALSE, cluster_cols = FALSE, 
                display_numbers = TRUE, number_format = "%.0f", 
                fontsize_number = 8)
        )
    }
}


#' @description \code{default.config.batch()} returns the path of 
#' basic configuration files for batch job submissions using \code{batchtools}.
#' These are stored in the internal package directory and allow working
#' with an SGE environment.
#' 
#' @param config.file character, of the form \code{batchtools.ist-X.Y}, where
#' \code{X} is \code{sge} or \code{email} and 
#' \code{Y} is \code{R} or \code{tmpl}
#' 
#' @return \code{default.config.batch()} returns a character with 
#' the path to the desired configuration file
#' 
#' @examples 
#' IST::default.config.batch("batchtools.ist-sge.tmpl")
#' 
#' @name defaults
#' @rdname defaults
#' 
#' @importFrom checkmate qassert test_subset
#' @export
default.config.batch <- function(config.file) {
    checkmate::qassert(config.file, "S1")
    checkmate::test_subset(
        config.file, 
        c(  "batchtools.ist-sge.R", 
            "batchtools.ist-email.R", 
            "batchtools.ist-sge.tmpl", 
            "batchtools.ist-email.tmpl"))
    
    system.file(paste0("templates/", config.file), package = "IST")
}
