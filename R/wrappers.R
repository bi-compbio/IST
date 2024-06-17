#' @title Wrapper functions
#' 
#' @description \code{do.isx()} runs the ISA/ISD/ISS flavour 
#' 
#' @param signature,mapping.data.table gene signature and (optional) 
#' precomputed orthology mapping
#' @param object.class character, either \code{"ist.discriminator"} to run the 
#' classical flavours (isa/iss/isd) or \code{"ist.pathways"} to fit
#' pathway models
#' @param name character, for traceability, will be returned as-is
#' @param org.from,org.to,X,y,de.genes,id.oc,id.bin,id.ref,id.ist arguments 
#' for in silico treatment
#' @param pathways.table.oc,pathways.table.bin arguments for pathway models
#' @param ... \code{do.isx()}: further arguments for fitting the discriminator, 
#' including the flavour of the ist; \code{do.batch.isx()}: ignored
#' 
#' @return \code{do.isx()}: \code{list} with the following fields: 
#' \code{name} for traceability, and \code{decisions}
#' 
#' @name wrappers
#' @rdname wrappers
#' 
#' @examples 
#' data(sample.data.ist)
#' ## IST flavours
#' list.flavours <- with(sample.data.ist, {
#' do.isx(
#' signature = fc.mmu.true, 
#' object.class = "ist.discriminator", 
#' mapping.data.table = df.orth, 
#' name = "test",
#' flavour = c("isd", "isa", "iss"),
#' org.from = "mmusculus", 
#' org.to = "hsapiens", 
#' X = X.hsa, 
#' y = y.hsa, 
#' id.oc = rep(TRUE, nrow(X.hsa)), 
#' id.bin = rep(TRUE, nrow(X.hsa)), 
#' id.ref = y.hsa == 1, 
#' id.ist = y.hsa == -1, 
#' ncomp = 3
#' )
#' })
#' summary(list.flavours$decisions)
#' ## IST pathways
#' list.pathways <- with(sample.data.ist, {
#' do.isx(
#' signature = fc.mmu.true, 
#' object.class = "ist.pathways", 
#' mapping.data.table = df.orth, 
#' name = "test",
#' pathways.table.oc = head(path.hsa.oc, 100),
#' pathways.table.bin = head(path.hsa.bin, 100),
#' org.from = "mmusculus", 
#' org.to = "hsapiens", 
#' X = X.hsa, 
#' y = y.hsa, 
#' id.oc = rep(TRUE, nrow(X.hsa)), 
#' id.bin = rep(TRUE, nrow(X.hsa)), 
#' id.ref = y.hsa == 1, 
#' id.ist = y.hsa == -1, 
#' ncomp = 3
#' )
#' })
#' summary(list.pathways$decisions)
#' 
#' @importFrom checkmate qassert assertSubset assertNames
#' @importFrom methods new
#' @importFrom plyr ldply 
#' @export do.isx
do.isx <- function(
    signature, object.class, 
    name, org.from, org.to, 
    X, y, ...) {
    
    checkmate::qassert(name, "S1")
    checkmate::qassert(object.class, "S1")
    checkmate::assertSubset(
        object.class, c("ist.discriminator", "ist.pathways"))
    
    args <- list(...)
    id.ref <- args$id.ref
    id.ist <- args$id.ist
    checkmate::assertNames(names(args), must.include = c("id.ref", "id.ist"))
    
    # instance and fit translator
    ist.trans <- define(
        "ist.translator", 
        animal.data.table = signature, 
        org.from = org.from, 
        org.to = org.to, 
        ...
    )
    ist.trans <- fit(ist.trans)
    
    # instance and fit discriminator
    ist.discr <- define(
        object.class, 
        X = X, 
        y = y, 
        org.to = org.to, 
        ...
    ) 
    ist.discr <- fit(ist.discr, ist.translator = ist.trans, ...)
    
    # transform human data (samples indicated in id.ist) with translator
    X.new <- predict(ist.trans, newdata = X[id.ist, , drop = FALSE])
    
    # manual prediction
    # prediction: pls scores isd/isa
    df.decisions <- plyr::ldply(
        list(trt = X.new, untrt = X),
        IST::predict, 
        object = ist.discr, 
        type = "decision",
        .id = "treatment")
    
    # keep track of metadata in data.frames as well
    df.decisions$sample.ist <- df.decisions$sample %in% 
        rownames(ist.discr@X)[id.ist]
    df.decisions$sample.ref <- df.decisions$sample %in% 
        rownames(ist.discr@X)[id.ref]
    df.decisions$name <- name
    
    list(
        name = name, 
        decisions = df.decisions
    )
}


#' @description \code{do.batch.isx()} submits an array of jobs that iterate
#' over `do.isx()`.
#' If \code{object.class == "ist.discriminator"}, then \code{de.genes} 
#' can be optionally set.
#' If \code{object.class == "ist.pathways"}, one should provide 
#' \code{pathways.table.bin} and/or \code{pathways.table.bin}
#' 
#' @param se.orgto \code{SummarizedExperiment} object with the \code{org.to}
#' molecular readouts (includes the variables \code{X} and \code{y}
#' implicitly)
#' @param list.sig list with all the signatures to apply in batch; 
#' names will be preserved in output
#' @param flavour character vector, flavours of ist
#' @param se.coldata.col character, column from the \code{SummarizedExperiment}
#' \code{colData()} that contains the response \code{y}
#' @param reg.name character, path to save the registry
#' @param waitForJobs logical, should the function wait until all the 
#' jobs are processed (and then apply a reduce)? If \code{FALSE}, 
#' the function returns right after the job submission.
#' @param conf.sge character, path to the \code{.R} configuration file for 
#' SGE (see \code{?makeClusterFunctionsSGE}) in batchtools
#' @param conf.email character, optional path to the \code{.R} 
#' configuration file for sending an e-mail after the main job is completed
#'  
#' @return \code{do.batch.isx()}: if \code{waitForJobs = TRUE}, it returns 
#' the reduced results, i.e. a data.frame with all the runs. Otherwise, 
#' it returns the \code{batchtools::batchMap} result
#' 
#' @name wrappers
#' @rdname wrappers
#' 
#' @importFrom SummarizedExperiment assay colData
#' @importFrom batchtools makeRegistry clearRegistry batchMap submitJobs
#' waitForJobs reduceResults sweepRegistry
#' @importFrom checkmate qassert
#' 
#' @export do.batch.isx
do.batch.isx <- function(
    se.orgto, list.sig, object.class, flavour,
    mapping.data.table, de.genes, 
    org.from, org.to, id.oc, id.bin, id.ref, id.ist, 
    pathways.table.oc, pathways.table.bin, 
    se.coldata.col = "ist.class", 
    reg.name = "reg-isx", 
    waitForJobs = TRUE,
    conf.sge = IST::default.config.batch("batchtools.ist-sge.R"), 
    conf.email, 
    ...) {
    
    checkmate::qassert(se.coldata.col, "S1")
    checkmate::qassert(reg.name, "S1")
    checkmate::qassert(waitForJobs, "B1")
    checkmate::qassert(conf.sge, "S1")
    
    # reference dataset
    X <- t(SummarizedExperiment::assay(se.orgto))
    y <- SummarizedExperiment::colData(se.orgto)[[se.coldata.col]]
    
    # set up batch job
    reg <- batchtools::makeRegistry(reg.name, conf.file = conf.sge)
    
    # make sure to clean registry
    batchtools::clearRegistry(reg)
    
    # do we have a registry?
    # getDefaultRegistry()
    
    # arguments
    # args <- list(signature = list.sig)
    more.args <- list(
        object.class = object.class, 
        org.from = org.from, org.to = org.to, 
        X = X, y = y, id.oc = id.oc, id.bin = id.bin, 
        id.ref = id.ref, id.ist = id.ist)
    
    # use de genes if avail
    if (!missing(de.genes)) 
        more.args$de.genes <- de.genes
    
    # use precomputed orthology mapping if available
    if (!missing(mapping.data.table)) 
        more.args$mapping.data.table <- mapping.data.table
    
    # add flavour or pathways.table.oc/bin if pathway models
    # TODO all thsese cases might work through ellipsis or list.args argument
    if (!missing(flavour))
        more.args$flavour <- flavour
    if (!missing(pathways.table.oc))
        more.args$pathways.table.oc <- pathways.table.oc
    if (!missing(pathways.table.bin))
        more.args$pathways.table.bin <- pathways.table.bin
    
    # define jobs
    jobMap <- batchtools::batchMap(
        fun = do.isx, 
        signature = list.sig, 
        name = names(list.sig), 
        more.args = more.args, 
        reg = reg)
    
    # see how they are not submitted yet
    # getStatus()
    
    # submit them
    batchtools::submitJobs(reg = reg)
    
    # This sys.sleep is not enough to assure some jobs have finished
    # Sys.sleep(15)
    
    # this breaks if the table is empty
    # View(findDone())
    
    # getStatus()
    
    # don't exit function until the job is completed
    if (waitForJobs) {
        batchtools::waitForJobs(reg = reg)
    } else {
        return(jobMap)
    }
    
    # send e-mail when it's over
    if (!missing(conf.email)) {
        checkmate::qassert(conf.email, "S1")
        
        dir.email <- paste0(tempdir(), "/email")
        # make sure no prior jobs left a registry
        unlink(dir.email, recursive = TRUE)
        
        reg.email <- batchtools::makeRegistry(
            dir.email, conf.file = conf.email)
        
        # dummy e-mail job
        batchtools::batchMap(fun = identity, args = list(0), reg = reg.email)
        batchtools::submitJobs(
            reg = reg.email, resources = list(job.name = "email-test"))
        batchtools::waitForJobs(reg = reg.email)
        batchtools::removeRegistry(wait = 0, reg = reg.email)
    }
    
    
    # also breaks if no errors are found
    # View(findErrors())
    # View(getJobStatus(reg = reg))
    
    # logs
    # tail(getLog(id = 3, reg = reg))
    # grepLogs(pattern = "Skipping", ignore.case = TRUE, reg = reg)
    
    # errors
    # getErrorMessages(reg = reg)
    
    # results of one job
    # lapply(1:5, function(x) head(loadResult(x, reg = reg)$decisions))
    
    # how to combine the results? 
    # bind the rows of the already combined i-1 jobs to the decisions of 
    # the i-th job
    f <- function(aggr, res) rbind(aggr, res$decisions)
    
    # make sure init is list(), otherwise it will apply 
    # f(job1, job2) on the first call, and we don't want that
    # We only want calls of the form
    # f(combined, jobk)
    # So I guess the first call is something like 
    # f(NULL, res) which works for us
    df.decisions <- batchtools::reduceResults(
        fun = f, init = list(), reg = reg)
    # dim(df.decisions)
    # rbind(head(df.decisions), tail(df.decisions))
    
    # removeRegistry(wait = 0, reg = reg)
    
    # read registry
    # regload <- loadRegistry("test-isa-batch")
    
    # Check Consistency and Remove Obsolete Information
    batchtools::sweepRegistry(reg)
    
    df.decisions
}
