context("Batch jobs for SGE")

skip_if_not_in_boehringer <- function() {
    # check if hostname looks like internal
    # logic taken out for publication
    if (TRUE) {
        skip("Only test batch functions inside Boehringer's servers")
    }
}

test_that("Batch IST", {
    skip("These batch functions are legacy material. Skipping tests...")
    
    skip_if_not_in_boehringer()
    
    data("sample.data.ist")
    set.seed(1)
    
    reg.name <- "reg-test-isx"
    if (dir.exists(reg.name)) 
        stop("Testing should not proceed if directory ", reg.name, " exists")
    
    org.from <- "mmusculus"
    org.to <- "hsapiens"
    id.oc <- rep(TRUE, nrow(sample.data.ist$X))
    id.bin <- rep(TRUE, nrow(sample.data.ist$X))
    id.ref <- sample.data.ist$y == 1
    id.ist <- sample.data.ist$y == -1
    mapping.data.table <- sample.data.ist$df.orth
    
    se.orgto <- SummarizedExperiment::SummarizedExperiment(
        assays = t(sample.data.ist$X.hsa), 
        colData = data.frame(ist.class = sample.data.ist$y.hsa)
    )
    # list with signatures
    list.sig <- sample.data.ist[c("fc.mmu.rnd", "fc.mmu.true")]
    de.genes <- head(colnames(sample.data.ist$X.hsa), 100)
    
    pathways.table.oc <- sample.data.ist$path.hsa.oc
    pathways.table.bin <- sample.data.ist$path.hsa.bin
    
    # now with email and predefined ortholog mapping
    expect_error({
        df.decisions <- do.batch.isx(
            se.orgto = se.orgto, 
            list.sig = list.sig, 
            object.class = "ist.discriminator", 
            flavour = c("isa", "isd", "iss"), 
            mapping.data.table = mapping.data.table,
            org.from = org.from, 
            org.to = org.to, 
            id.oc = id.oc, id.bin = id.bin, 
            id.ref = id.ref, id.ist = id.ist, 
            # de.genes = de.genes,
            se.coldata.col = "ist.class", 
            reg.name = reg.name)
            # , 
            # conf.email = default.config.batch("batchtools.ist-email.R")
            # )
    }, NA)
    unlink(reg.name, recursive = TRUE)
    expect_is(df.decisions, "data.frame")
    
    # make sure all flavours were run
    expect_setequal(unique(df.decisions$flavour), c("isa", "isd", "iss"))
    
    # try pathway models
    expect_error({
        df.decisions <- do.batch.isx(
            se.orgto = se.orgto, 
            list.sig = list.sig, 
            object.class = "ist.pathways", 
            pathways.table.oc = pathways.table.oc, 
            pathways.table.bin = pathways.table.bin, 
            mapping.data.table = mapping.data.table,
            org.from = org.from, 
            org.to = org.to, 
            id.oc = id.oc, id.bin = id.bin, 
            id.ref = id.ref, id.ist = id.ist, 
            se.coldata.col = "ist.class", 
            reg.name = reg.name)
    }, NA)
    unlink(reg.name, recursive = TRUE)
    expect_is(df.decisions, "data.frame")
    
    # make sure all the pathways in the input were run
    expect_setequal(
        unique(pathways.table.oc$path.id), 
        unique(df.decisions$path.id))
    
    # Info 
    # 
    # SGE does not like when the registry is in a temporary location
    # (cannot find log file, even when the directory is correct)
    # 
    # to kill all jobs: 
    # # batchtools::killJobs(batchtools::getJobTable()$job.id)
    # to remove registry: 
    # # batchtools::removeRegistry(wait = 0)
    # 
    # Terminal: 
    # to follow status: watch qstat
    # to delete all my jobs: qdel -u $(whoami)
    # to get error msgs: qstat -j JOB_ID | grep error
    # 
    # To debug
    # # load results
    # reg <- batchtools::loadRegistry(reg.name, writeable = TRUE)
    # # get status
    # reg$status
    # # remove registry and move on
    # removeRegistry(wait = 0)
})
