library(IST)
data("sample.data.ist")

X <- sample.data.ist$X.hsa
y <- sample.data.ist$y.hsa
ref.y <- 1

sigrnd <- sample.data.ist$fc.mmu.rnd
sigtrue <- sample.data.ist$fc.mmu.true

org.from <- "mmusculus"
org.to <- "hsapiens"

do.isa <- function(X, y, ref.y, org.from, signature, org.to, ...) {
    # instance and fit translator
    ist.trans <- IST:::ist.translator(
        animal.data.table = signature, 
        org.from = org.from, 
        org.to = org.to) 
    
    ist.trans <- fit(ist.trans)
    
    # dummy line: need gene names, even if not used
    de.genes <- head(colnames(X), 100)
    
    # transform human data
    X.new <- predict(ist.trans, newdata = X)
    
    # instance and fit discriminator
    ist.discr <- IST:::ist.discriminator(
        X = X, 
        y = y, 
        org.to = org.to, 
        de.genes = de.genes)
    
    ist.discr <- fit(
        ist.discr, 
        flavour = c("isa"), 
        ref.y = ref.y, 
        ist.translator = ist.trans, 
        ncomp = 10)
    
    # manual prediction
    # prediction: pls scores isd/isa
    df.scores <- plyr::ldply(
        list(trt = X.new, untrt = X),
        predict, 
        object = ist.discr, 
        type = "scores",
        .id = "treatment")
    df.decisions <- plyr::ldply(
        list(trt = X.new, untrt = X),
        predict, 
        object = ist.discr, 
        type = "decision",
        .id = "treatment")
    
    list(scores = df.scores, 
         decisions = df.decisions, 
         isa = ist.discr@isa)
}

res.isa <- do.isa(
    X = X, y = y, ref.y = ref.y, org.from = org.from, 
    signature = sigtrue, org.to = org.to
)

ggplot(res.isa$scores, 
       aes(x = Comp1, y = Comp2, colour = gsub("\\d", "", sample))) +
    geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
    geom_vline(xintercept = 0, lty = 2, colour = "gray50") +
    geom_point() +
    facet_grid(flavour~treatment) + 
    theme_bw() +
    theme(aspect.ratio = 1)

ggplot(dplyr::filter(res.isa$scores, treatment == "trt"), 
       aes(x = Comp1, y = Comp2, colour = gsub("\\d", "", sample))) +
    geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
    geom_vline(xintercept = 0, lty = 2, colour = "gray50") +
    stat_density_2d(
        aes(x = Comp1, y = Comp2, colour = gsub("\\d", "", sample), 
            fill = stat(level)), geom = "polygon", 
        data = dplyr::filter(res.isa$scores, treatment == "untrt")) +
    scale_fill_gradient(low = "gray90", high = "black") +
    geom_point(size = 1, pch = 19) +
    facet_grid(~flavour) + 
    theme_bw() +
    theme(aspect.ratio = 1)

ggplot(res.isa$decisions, aes(x = gsub("\\d", "", sample), y = decision.value)) +
    geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
    geom_boxplot() +
    facet_grid(flavour~treatment) + 
    theme_bw() +
    theme(aspect.ratio = 1)
