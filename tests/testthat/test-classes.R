context("Classes")

test_that("Translator", {
    expect_error({
        ist.translator(
            animal.data.table = iris[, 1:3], 
            org.from = "a", org.to = "b")
    }, "names")
    
    dat <- cars
    colnames(dat) <- c("logFC" ,"adj.P.Val")
    
    expect_error({
        ist.translator(
            animal.data.table = dat, 
            org.from = "a", org.to = "b"
        )
    }, NA)
    
    # rownames(dat) <- NULL
    # expect_error({
    #     ist.translator(
    #         animal.data.table = dat, 
    #         org.from = "a", org.to = "b"
    #     )
    # }, "rownames")

})

test_that("Signatures", {
    data(sample.env.ist)
    
    expect_error({
        ist.sig <- with(
            sample.env.ist, 
            define(
                "ist.signatures", 
                list.sig = list.sig, 
                tab.meta = df.meta.sig, 
                org.to = "hsapiens"
            )
        )
    }, NA)
})

test_that("Discriminator", {
    X <- matrix(1:120, nrow = 3)
    y <- c(0, 0, 1)
    id.oc <- id.bin <- id.ref <- id.ist <- rep(TRUE, nrow(X))
    
    # X must be row/colnamed
    expect_error({
        define(
            "ist.discriminator", 
            X = X, y = y, org.to = "b"
        )
    }, "names")
    
    rownames(X) <- seq_len(nrow(X))
    colnames(X) <- seq_len(ncol(X))
    
    # correct format
    expect_error({
        define(
            "ist.discriminator", 
            X = X, y = y, org.to = "b", 
            id.oc = id.oc, id.bin = id.bin, 
            id.ref = id.ref, id.ist = id.ist, 
            de.genes = as.character(1:10)
        )
    }, NA)
    
    # incorrect format for de.genes
    expect_error({
        define(
            "ist.discriminator", 
            X = X, y = y, org.to = "b", de.genes = 1:10
        )
    }, "de.genes")
    
    # mismatching dimensions
    expect_error({
        define(
            "ist.discriminator", 
            X = X, y = 1, org.to = "b"
        )
    }, "as many rows")
    
    # missings in X
    expect_error({
        define(
            "ist.discriminator", 
            X = cbind(X, NAs = NA), y = 1, org.to = "b"
        )
    }, "missing")
    
    # missings in y
    expect_error({
        define(
            "ist.discriminator", 
            X = X, y = c(1, 1, NA), org.to = "b"
        )
    }, "missing")
})

test_that("Pathways", {
    X <- matrix(1:120, nrow = 3)
    y <- c(0, 0, 1)
    id.oc <- id.bin <- id.ref <- id.ist <- rep(TRUE, nrow(X))
    
    path.oc <- data.frame(path = as.character(1:10), gene = as.character(1:10), 
                          stringsAsFactors = FALSE)
    path.bin <- path.oc
    
    # X must be row/colnamed
    expect_error({
        define(
            "ist.pathways", 
            X = X, y = y, org.to = "b"
        )
    }, "names")
    
    rownames(X) <- seq_len(nrow(X))
    colnames(X) <- seq_len(ncol(X))
    
    # correct format
    define(
        "ist.pathways", 
        X = X, y = y, org.to = "b", id.oc = id.oc, id.bin = id.bin, 
        pathways.table.oc = path.oc, pathways.table.bin = path.bin
    )
})
