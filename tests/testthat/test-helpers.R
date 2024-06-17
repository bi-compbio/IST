context("Helper functions")

test_that("basic limma", {
    data(sample.data.ist)
    
    # numeric
    # warning from limma
    expect_warning({
        de.genes.num <- get.de.genes(
            sample.data.ist$X.hsa, 
            sample.data.ist$y.hsa, 
            return.toptable = TRUE
        )
        
    }, "row names")
    
    checkmate::expect_data_frame(de.genes.num$topTable)
    checkmate::expect_character(de.genes.num$de.genes)
    
    # integer
    # warning from limma
    expect_warning({
        de.genes.int <- get.de.genes(
            sample.data.ist$X.hsa, 
            as.integer(sample.data.ist$y.hsa), 
            return.toptable = TRUE
        )
    }, "row names")
    
    checkmate::expect_data_frame(de.genes.int$topTable)
    checkmate::expect_character(de.genes.int$de.genes)
    
    # factor
    # warning from limma
    expect_warning({
        de.genes.fac <- get.de.genes(
            sample.data.ist$X.hsa, 
            factor(sample.data.ist$y.hsa, levels = c(-1, 1)), 
            return.toptable = TRUE
        )
    }, "row names")
    
    checkmate::expect_data_frame(de.genes.fac$topTable)
    checkmate::expect_character(de.genes.fac$de.genes)
    
    # numeric and int should be the same thing
    expect_equivalent(de.genes.num, de.genes.int)
    
    # numeric and factor should have a scaling constant of 1/2
    # (numeric -1/1, factor 0/1)
    expect_equivalent(
        de.genes.num$topTable$logFC, 
        de.genes.fac$topTable$logFC/2)
    
    # ..and the same DE genes
    expect_equivalent(de.genes.num$de.genes, de.genes.fac$de.genes)
})
