context("Default arguments")

test_that("default arguments match method formals", {
    # which arguments are required for train.default?
    tr.method <- getAnywhere("train.default")   
    tr.formals <- names(formals(tr.method$objs[[1]]))
    
    expect_error({
        arg.train <- default.args.bin.train()
    }, NA)
    
    expect_true(
        all(names(arg.train) %in% tr.formals)
    )
    
    # same with pls
    pls.formals <- c(
        names(formals(pls::mvr)), 
        names(formals(pls::crossval))
    )
    
    expect_error({
        arg.pls <- default.args.bin.pls()
    }, NA)
    
    expect_true(
        all(names(arg.pls) %in% pls.formals)
    )
})
