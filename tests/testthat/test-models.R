context("Models")

test_that("One-class model", {
    data("sample.data.ist")
    
    set.seed(1)
    X.de <- head(
        sample(colnames(sample.data.ist$X.hsa)), 
        50
    )
    
    expect_error({
        set.seed(123)
        mod.iss <- fit.ist.oc(
            X.ref = sample.data.ist$X.hsa, 
            de.genes = X.de, 
            n.genes = 100
        )
        
        set.seed(234)
        mod.iss2 <- fit.ist.oc(
            X.ref = sample.data.ist$X.hsa, 
            de.genes = X.de, 
            n.genes = 100
        )
    }, NA)
    
    # seed is set internally, so external seed should not change anything
    expect_equivalent(mod.iss, mod.iss2)
    
    expect_error({
        pred.iss <- predict.ist.oc(
            mod.iss, 
            newdata = sample.data.ist$X.hsa
        )
    }, NA)
    
    expect_error({
        pred.iss <- predict.ist.oc(mod.iss)
    }, "newdata")

})

test_that("Binary classification model", {
    data("sample.data.ist")
    
    set.seed(1)
    X.de <- head(
        sample(colnames(sample.data.ist$X.hsa)), 
        50
    )
    
    expect_error({
        set.seed(2)
        mod.pls <- fit.ist.bin(
            X = sample.data.ist$X.hsa, 
            y = sample.data.ist$y.hsa,
            de.genes = X.de
        )
        
        set.seed(3)
        mod.pls2 <- fit.ist.bin(
            X = sample.data.ist$X.hsa, 
            y = sample.data.ist$y.hsa,
            de.genes = X.de
        )
    }, NA)
    
    # has the data been centered and scaled?
    expect_true({
        all(c("center", "scale") %in% names(mod.pls@train$preProcess$method))
    })
    
    # expect no differences between runs exept "times" (execution times)
    # 
    # fields that should not change
    nm.train <- setdiff(names(mod.pls@train), c("times", "finalModel"))
    
    expect_equivalent(mod.pls@train[nm.train], mod.pls2@train[nm.train])
    
    # must store optimal Ncomp
    expect_equal(
        mod.pls@train$finalModel$ncomp, 
        mod.pls@Ncomp
    )
    
    # in this synthetic example, the optimal number of components (around 2) 
    # is less than the total amount of fitted ones (5)
    expect_gt(
        mod.pls@pls.mod$ncomp, 
        mod.pls@Ncomp
    )
    
    expect_error({
        pred.pls <- predict.ist.bin(
            mod.pls, 
            newdata = sample.data.ist$X.hsa
        )
    }, NA)
    
    expect_error({
        pred.pls <- predict.ist.bin(mod.pls)
    }, "newdata")
    
})
