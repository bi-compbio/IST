default.resources <- list(walltime = 3600,
                          queue = "tdr.rh7",
                          measure.memory = TRUE)

# max.concurrent.jobs <- 2

# Now just try to provide the path using IST's helper!
cluster.functions <- batchtools::makeClusterFunctionsSGE(
    template = IST::default.config.batch("batchtools.ist-email.tmpl"))
