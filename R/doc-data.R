#' Sample data to run IST
#' 
#' The sample data is a list with the 
#' essential elements to run the complete IST.
#' Comes from synthetic gaussian random variables
#' 
#' @format A list with human expression \code{X.hsa}, 
#' human response \code{y.hsa}, 
#' murine random fold changes \code{fc.mmu.rnd}, 
#' (some of the) murine true fold changes \code{fc.mmu.true}, 
#' \code{path.hsa.oc} and \code{path.hsa.bin} for synthetic pathways 
#' (random fixed-size ones) and, 
#' for traceability, \code{df.orth} with 
#' the orthology mappings within 
#' the sample data.
#' For the multi-signature analysis, 
#' the slot \code{dt.fc.mmu} contains murine signatures and 
#' \code{df.meta.mmu} their metadata
"sample.data.ist"

#' Sample data to run IST (II)
#' 
#' The sample data is an environment with the 
#' essential elements to run the complete IST.
#' This is an updated version of \code{sample.data.ist}, also using
#' synthetic gaussian random variables.
#' 
#' @format An environment with human expression \code{X.hsa}, 
#' human response \code{y.hsa}, list of signatures \code{list.sig} 
#' and their metadata \code{df.meta.sig}, table of pathways 
#' \code{df.path} and their metadata \code{df.meta.path}, 
#' and other non-essential variables (\code{de.genes} for 
#' a list of differential limma genes and \code{list.orth} for 
#' the orthology mappings)
"sample.env.ist"

#' Sample results object from IST run
#' 
#' The sample output data is an \code{ist.results} object, 
#' allows playing with plotting and interactive visualisation.
#' Comes from the synthetic data in \code{sample.env.ist}
#' 
#' @format An object of class \code{ist.results} which implicitly
#' contains the \code{ist.pathways} and the \code{ist.signatures} objects
#' used to fit it. 
"sample.results.ist"

#' Default orthology mapping from mouse/rat/human to human
#' 
#' Orthology mappings that can be plugged into 
#' \code{ist.translator} and \code{ist.signatures} objects.
#' These were generated using \code{orthIDcon}.
#' Generated March 30th, 2020 using biomaRt.
#' 
#' @format A \code{list} with names hsapiens, mmusculus and rnorvegicus, 
#' with \code{data.table}s as provided by \code{orthIDcon}
"data.list.orth"

#' Default labels for human ENSEMBL genes
#' 
#' These mappings (\code{vec.ensembl2symbol}, \code{vec.ensembl2symbolonly},
#' \code{vec.ensembl2entrez}) allow a convenient mapping
#' between ENSEMBL identifiers and more understandable labels
#' when plotting an \code{ist.results} object.
#' Generated March 30th, 2020 using biomaRt.
#' 
#' @format A named vector, directly as needed by \code{ist.results}
#' 
#' @aliases vec.ensembl2symbol vec.ensembl2symbolonly vec.ensembl2entrez
"vec.ensembl2symbol"

