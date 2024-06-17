#' @title Fitting methods
#' 
#' @include fit.R
#' @include predict.R
#' @include show.R
#' @include plot.R
#' 
#' @name Fit
#' @rdname fit-methods
methods::setGeneric(
    "fit", 
    def = function(object, ...) {})

#' @importFrom methods signature
#' @exportMethod fit
methods::setMethod(
    "fit", 
    signature = methods::signature(object = "ist.translator"), 
    fit.ist.translator)

#' @importFrom methods signature
#' @exportMethod fit
methods::setMethod(
    "fit", 
    signature = methods::signature(object = "ist.signatures"), 
    fit.ist.signatures)

#' @importFrom methods signature
#' @exportMethod fit
methods::setMethod(
    "fit", 
    signature = methods::signature(object = "ist.discriminator"), 
    fit.ist.discriminator)

#' @importFrom methods signature
#' @exportMethod fit
methods::setMethod(
    "fit", 
    signature = methods::signature(object = "ist.pathways"), 
    fit.ist.pathways)

#' @importFrom methods signature
#' @exportMethod fit
methods::setMethod(
    "fit", 
    signature = methods::signature(object = "ist.results"), 
    fit.ist.results)




#' @title Predict methods
#' 
#' @name Predict
#' @rdname predict-methods
methods::setGeneric("predict")

#' @importFrom methods signature
#' @exportMethod predict
methods::setMethod(
    "predict", 
    signature = methods::signature(object = "ist.translator"), 
    predict.ist.translator)

#' @importFrom methods signature
#' @exportMethod predict
methods::setMethod(
    "predict", 
    signature = methods::signature(object = "ist.signatures"), 
    predict.ist.signatures)

#' @importFrom methods signature
#' @exportMethod predict
methods::setMethod(
    "predict", 
    signature = methods::signature(object = "ist.mod.oc"), 
    predict.ist.oc)

#' @importFrom methods signature
#' @exportMethod predict
methods::setMethod(
    "predict", 
    signature = methods::signature(object = "ist.mod.bin"), 
    predict.ist.bin)

#' @importFrom methods signature
#' @exportMethod predict
methods::setMethod(
    "predict", 
    signature = methods::signature(object = "ist.discriminator"), 
    predict.ist.discriminator)

#' @importFrom methods signature
#' @exportMethod predict
methods::setMethod(
    "predict", 
    signature = methods::signature(object = "ist.pathways"), 
    predict.ist.pathways)





#' @title Show methods
#'
#' @name Show
#' @importFrom methods signature
#' @exportMethod show
methods::setMethod(
    "show", 
    signature = methods::signature(object = "ist.translator"), 
    show.ist.translator)

#' @importFrom methods signature
#' @exportMethod show
methods::setMethod(
    "show", 
    signature = methods::signature(object = "ist.signatures"), 
    show.ist.signatures)

#' @importFrom methods signature
#' @exportMethod show
methods::setMethod(
    "show", 
    signature = methods::signature(object = "ist.discriminator"), 
    show.ist.discriminator)

#' @importFrom methods signature
#' @exportMethod show
methods::setMethod(
    "show", 
    signature = methods::signature(object = "ist.pathways"), 
    show.ist.pathways)

#' @importFrom methods signature
#' @exportMethod show
methods::setMethod(
    "show", 
    signature = methods::signature(object = "ist.results"), 
    show.ist.results)



#' @title Plot methods
#' 
#' @name Plot
#' @importFrom methods signature
#' @exportMethod plot
methods::setMethod(
    "plot", 
    signature = methods::signature(
        x = "ist.translator", 
        y = "ist.discriminator"), 
    plot.ist.all)
