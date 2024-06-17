default.resources <- list(walltime = 3600,
                          queue = "tdr.rh7",
                          measure.memory = TRUE)

# max.concurrent.jobs <- 2

dir.tmpl <- IST::default.config.batch("batchtools.ist-sge.tmpl")

message("Using the following configuration file: ", dir.tmpl)

# The template is given to be found in the following way
# (1st option in ?makeClusterFunctionsSGE)
# 
# “batchtools.[template].tmpl” in the path specified by 
# the environment variable “R_BATCHTOOLS_SEARCH_PATH”.
# 
# ... no, instead just try to provide the path using IST's helper!
cluster.functions <- batchtools::makeClusterFunctionsSGE(template = dir.tmpl)
