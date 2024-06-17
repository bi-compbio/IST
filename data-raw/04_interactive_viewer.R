library(shiny)

# set up paths for rds object with ist.results, app title and app description
path_current <- getwd()
shinyOptions(ist.result.path = paste0(path_current, "/data-raw/01_toy_example/03_data/sample.results.ist.rds"))
shinyOptions(ist.browser.title = "Sample data")
shinyOptions(ist.browser.description = "Sample data. 30 patients, 30 controls.")

# need the ISTBrowser package to run the app (either installed, or loaded)
# devtools::load_all("../istbrowser")
shiny::shinyAppDir(system.file("app", package = "ISTBrowser"))
