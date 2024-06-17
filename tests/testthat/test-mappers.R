context("Mapping functions")

# to avoid proxy errors
skip_if_not_in_boehringer <- function() {
  # check if hostname looks like a BI machine
  if (!grepl("boehringer.com", Sys.info()["nodename"])) {
    skip(c(
      "Only test batch functions inside Boehringer's servers;", 
      "the bundled SGE config files are not expected to work outside."))
  }
}

test_that("Orthology mapping", {
    skip_if_not_in_boehringer()
  
    expect_error({
        df.human <- orthIDsingle("hsapiens")
    }, NA)
    
    expect_gt(nrow(df.human), 50000)
  
    expect_error({
        df.mouse2human <- orthIDcon("mmusculus", "hsapiens")
    }, NA)
  
    expect_gt(nrow(df.mouse2human), 10000)
})
