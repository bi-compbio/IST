context("Workflow")

test_that("Fitting translator and discriminator", {
    data("sample.data.ist")
    data("sample.env.ist")
    set.seed(1)
    
    de.genes <- head(sample(colnames(sample.data.ist$X.hsa)), 100)
    
    # fit translator with predefined mapping
    expect_error({
        ist.trans.pre <- ist.translator(
            animal.data.table = sample.data.ist$fc.mmu.rnd, 
            mapping.data.table = sample.data.ist$df.orth,
            org.from = "mmusculus", 
            org.to = "hsapiens")
    }, NA)
    
    expect_message({
        ist.trans.pre <- fit(ist.trans.pre)
    }, "Mapped")
    
    expect_error({
        X.new.pre <- predict(ist.trans.pre, newdata = sample.data.ist$X.hsa)
    }, NA)
    
    # fit translator without predefined mapping
    expect_error({
        ist.trans <- ist.translator(
            animal.data.table = sample.data.ist$fc.mmu.rnd, 
            org.from = "mmusculus", 
            org.to = "hsapiens")
    }, NA)
    
    expect_message({
        ist.trans <- fit(ist.trans)
    }, "Mapped")
    
    expect_error({
        X.new <- predict(ist.trans, newdata = sample.data.ist$X.hsa)
    }, NA)
    
    # fit multiple signatures
    expect_error({
        ist.signatures <- local({
            define(
                "ist.signatures", 
                list.sig = list.sig, tab.meta = df.meta.sig, 
                org.to = "hsapiens", list.mapping = list.orth
            )}, sample.env.ist
        )
        
        ist.signatures <- fit(ist.signatures)
    }, NA)
    
    # predict on multiple signatures
    expect_error({
        pred.signatures <- predict(
            ist.signatures, newdata = sample.env.ist$X.hsa)
    }, NA)
    expect_equal(names(pred.signatures), names(sample.env.ist$list.sig))
    
    # discriminator
    expect_error({
        ist.discr <- ist.discriminator(
            X = sample.data.ist$X, 
            y = sample.data.ist$y, 
            id.oc = rep(TRUE, nrow(sample.data.ist$X)), 
            id.bin = rep(TRUE, nrow(sample.data.ist$X)), 
            id.ref = sample.data.ist$y == 1, 
            id.ist = sample.data.ist$y == -1, 
            org.to = "hsapiens", 
            de.genes = de.genes)
    }, NA)
    
    expect_error({
        ist.discr <- fit(
            ist.discr, 
            flavour = c("iss", "isd", "isa"), 
            ist.translator = ist.trans, 
            ncomp = 10)
    }, NA)
    
    expect_error({
        decision.discr <- predict(ist.discr, newdata = X.new, type = "decision")
    }, NA)
    expect_true("decision" %in% colnames(decision.discr))
    
    expect_error({
        scores.discr <- predict(ist.discr, newdata = X.new, type = "scores")
    }, NA)
    expect_true("Comp1" %in% colnames(scores.discr))
    
    expect_error({
        plot(ist.trans, ist.discr, type = "decision")
        plot(ist.trans, ist.discr, type = "scores")
    }, NA)
    
    # discriminator that finds DE genes
    expect_error({
        ist.discr.nogenes <- ist.discriminator(
            X = sample.data.ist$X, 
            y = sample.data.ist$y, 
            id.oc = rep(TRUE, nrow(sample.data.ist$X)), 
            id.bin = rep(TRUE, nrow(sample.data.ist$X)), 
            id.ref = sample.data.ist$y == 1, 
            id.ist = sample.data.ist$y == -1, 
            org.to = "hsapiens")
    }, NA)
    
    # warning from limma
    expect_warning({
        ist.discr.nogenes <- fit(
            ist.discr.nogenes, 
            flavour = c("isa"), 
            ist.translator = ist.trans, 
            ncomp = 2)
    }, "row names")
    
    # Things that shold raise messages/warnings
    # 
    # Incomplete is* should warn and should predict correctly
    # but only on the present ist flavours
    expect_error({
        ist.discr.noiss <- ist.discriminator(
            X = sample.data.ist$X, 
            y = sample.data.ist$y, 
            id.oc = rep(TRUE, nrow(sample.data.ist$X)), 
            id.bin = rep(TRUE, nrow(sample.data.ist$X)), 
            id.ref = sample.data.ist$y == 1, 
            id.ist = sample.data.ist$y == -1, 
            org.to = "hsapiens", 
            de.genes = de.genes)
    
        ist.discr.noiss <- fit(
            ist.discr.noiss, 
            flavour = "isd", 
            ist.translator = ist.trans, 
            ncomp = 10)
    }, NA)
    
    expect_message({
        decision.discr.noiss <- predict(ist.discr.noiss, newdata = X.new)
    }, "iss")
})


test_that("Fitting signatures and pathways", {
    data("sample.env.ist")
    data("vec.ensembl2symbol")
    set.seed(1)
    
    ist.sig <- local({
        define(
            "ist.signatures", 
            list.sig = list.sig, 
            tab.meta = df.meta.sig,
            org.to = "hsapiens",
            list.mapping = list.orth
        )
    }, sample.env.ist)
    
    ist.sig <- fit(ist.sig)
    
    ist.path <- local({
        define(
            "ist.pathways", 
            X = X.hsa, 
            y = y.hsa, 
            org.to = "hsapiens", 
            id.oc = rep(TRUE, nrow(X.hsa)), 
            id.bin = rep(TRUE, nrow(X.hsa)), 
            pathways.table.oc = df.path, 
            pathways.table.bin = df.path, 
            pathways.meta = df.meta.path
        )
    }, sample.env.ist)
    
    # for some reason, executing fit() throws warning, but it is supressed
    # when wrapped in expect_warning
    suppressWarnings({
        ist.path <- fit(ist.path)
    })
    
    expect_s4_class(ist.path@list.oc[[1]], "ist.mod.oc")
    expect_s4_class(ist.path@list.bin[[1]], "ist.mod.bin")
    
    df.scores.path <- predict(
        ist.path, sample.env.ist$X.hsa, type = "scores"
    )
    
    df.decision.path <- predict(
        ist.path, sample.env.ist$X.hsa, type = "decision"
    )
    
    expect_s3_class(df.scores.path, "data.frame")
    expect_s3_class(df.decision.path, "data.frame")

    ist.res <- define(
        "ist.results", 
        ist.signatures = ist.sig, 
        ist.pathways = ist.path, 
        id.ref = sample.env.ist$y.hsa == 1, 
        id.ist = sample.env.ist$y.hsa == -1, 
        vec.gene2label = vec.ensembl2symbol 
    )
    
    ist.res <- fit(ist.res)
    
    path.nm <- getPathways(ist.res)
    
    # plots
    # 
    # boxplots
    plot.ist.boxplots(ist.res, y = path.nm, mapping = aes(fill = sig.org))
    
    plot.ist.boxplots(
        ist.res, y = path.nm, mapping = aes(fill = sig.org), 
        facet_rows = vars(sig.org), main = "Plot all the pathways!", 
        main.width = 15)
    
    # gene heatmaps
    plot.ist.genemaps(ist.res, y = "signal2", text.size = 2)
    
    plt.gene <- plot.ist.genemaps(
        ist.res, y = "signal1", type = "pheatmap", 
        main = "here, a really really really really really really long title",
        vars.meta.sig = c("sig.type", "sig.org"))
    
    # limit number of genes to 1 (extreme case)
    # check pheatmap
    plt.gene.trim.ph <- plot.ist.genemaps(
        ist.res, y = "signal1", type = "pheatmap", max.genes = 1,
        main = "here, a really really really really really really long title",
        vars.meta.sig = c("sig.type", "sig.org"))
    expect_equal(ncol(plt.gene.trim.ph$plot.data$data.wide), 1)
    # check ggplot
    plt.gene.trim.gg <- plot.ist.genemaps(
        ist.res, y = "signal1", type = "ggplot", max.genes = 1,
        main = "here, a really really really really really really long title",
        vars.meta.sig = c("sig.type", "sig.org"))
    expect_equal(length(unique(plt.gene.trim.gg$plot.data$ortholog)), 1)
    
    # pathway heatmaps
    plot.ist.pathwaymaps(ist.res, text.size = 5)
    # pathway heatmaps of top pathways and signatures only
    plot.ist.pathwaymaps(
        ist.res, 
        sig.ids = rankPathwaymaps(ist.res, what = "signatures", max.out = 3), 
        y = rankPathwaymaps(ist.res, what = "pathways", max.out = 3), 
        type = "pheatmap", 
        vars.meta.sig = c("sig.type", "sig.org"), 
        vars.meta.path = c("path.source", "path.name"), 
        main = "A good overall summary!")
    
    # save plots
    dir.out <- tempdir()
    file.out <- tempfile()
    save.boxplots(ist.res, paste0(dir.out, "/boxplots"))
    
    save.pheatmap(plt.gene, filename = paste0(file.out, "-p.png"))
    
    save.genemaps(ist.res, paste0(dir.out, "/genemaps"))
    # limit to 5 genes everywhere
    save.genemaps(ist.res, paste0(dir.out, "/genemaps-top5"), max.genes = 5)
    # browse to temp folder dir.out in rstudio to check
    
    # save.genemaps(ist.res, ".", sig.ids = head(getSignatures(ist.res), 1))
    # should not save anything
    save.genemaps(ist.res, paste0(dir.out, "/genemaps"), sig.ids = "inexistent")
})
