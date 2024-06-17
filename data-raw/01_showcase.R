devtools::load_all()
library(ggplot2)
library(pcaMethods)

library(magrittr)
library(plyr)
library(dplyr)

data("sample.data.ist")
set.seed(1)

df.pca <- pca(sample.data.ist$X.hsa, scale = "uv", center = TRUE) %>% 
    scores %>%
    as.data.frame %>%
    dplyr::mutate(sample = rownames(.), class = gsub("\\d", "", sample))

ggplot(df.pca, aes(x = PC1, y = PC2, colour = class)) +
    geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
    geom_vline(xintercept = 0, lty = 2, colour = "gray50") +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio = 1)

de.genes <- head(sample(colnames(sample.data.ist$X.hsa)), 100)

# p1: bad translator
ist.trans <- ist.translator(
    animal.data.table = sample.data.ist$fc.mmu.rnd, 
    org.from = "mmusculus", 
    org.to = "hsapiens")

ist.trans <- fit(ist.trans)

# p1: good translator
ist.trans <- ist.translator(
    animal.data.table = sample.data.ist$fc.mmu.true, 
    org.from = "mmusculus", 
    org.to = "hsapiens")

ist.trans <- fit(ist.trans)

X.new <- predict(ist.trans, newdata = sample.data.ist$X.hsa)

ist.discr <- ist.discriminator(
    X = sample.data.ist$X, 
    y = sample.data.ist$y, 
    org.to = "hsapiens", 
    de.genes = de.genes)

ist.discr <- fit(
    ist.discr, 
    flavour = c("iss", "isd", "isa"), 
    ref.y = 1, 
    ist.translator = ist.trans, 
    ncomp = 10)

# default method: prediction + plot
plot(ist.trans, ist.discr, type = "scores")
plot(ist.trans, ist.discr, type = "decision")

# pathway models
ist.path <- ist.pathways(
    X = sample.data.ist$X, 
    y = sample.data.ist$y, 
    org.to = "hsapiens", 
    ref.y = 1,
    pathways.table.oc = sample.data.ist$path.hsa.oc %>% head(100), 
    pathways.table.bin = sample.data.ist$path.hsa.bin %>% head(100)
)
ist.path <- fit(ist.path)

df.scores.path <- plyr::ldply(
    list(trt = X.new, untrt = sample.data.ist$X.hsa),
    predict, 
    object = ist.path, 
    type = "scores",
    .id = "treatment")

df.decision.path <- plyr::ldply(
    list(trt = X.new, untrt = sample.data.ist$X.hsa),
    predict, 
    object = ist.path, 
    type = "decision",
    .id = "treatment")


# manual prediction
# prediction: pls scores isd/isa
df.scores <- plyr::ldply(
    list(trt = X.new, untrt = sample.data.ist$X.hsa),
    predict, 
    object = ist.discr, 
    type = "scores",
    .id = "treatment")

ggplot(df.scores, aes(x = Comp1, 
                      y = Comp2, 
                      colour = gsub("\\d", "", sample))) +
    geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
    geom_vline(xintercept = 0, lty = 2, colour = "gray50") +
    geom_point() +
    facet_grid(flavour~treatment) + 
    theme_bw() +
    theme(aspect.ratio = 1)

ggplot(dplyr::filter(df.scores, treatment == "trt"), 
       aes(x = Comp1, 
           y = Comp2, 
           colour = gsub("\\d", "", sample))) +
    geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
    geom_vline(xintercept = 0, lty = 2, colour = "gray50") +
    stat_density_2d(aes(x = Comp1, 
                        y = Comp2, 
                        colour = gsub("\\d", "", sample), 
                        fill = stat(level)), geom= "polygon", 
                    data = dplyr::filter(df.scores, treatment == "untrt")) +
    scale_fill_gradient(low = "gray90", high = "black") +
    geom_point(size = 1, pch = 19) +
    facet_grid(~flavour) + 
    theme_bw() +
    theme(aspect.ratio = 1)

# prediction: decision values (all flavours)
df.decision <- plyr::ldply(
    list(trt = X.new, untrt = sample.data.ist$X.hsa),
    predict, 
    object = ist.discr, 
    .id = "treatment")

ggplot(df.decision, aes(x = gsub("\\d", "", sample), 
                        y = decision.value)) +
    geom_hline(yintercept = 0, lty = 2, colour = "gray50") +
    geom_boxplot() +
    facet_grid(flavour~treatment) + 
    theme_bw() +
    theme(aspect.ratio = 1)
