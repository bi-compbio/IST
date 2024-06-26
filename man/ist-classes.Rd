% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{ist.mod.oc-class}
\alias{ist.mod.oc-class}
\alias{ist.mod.oc}
\alias{ist.mod.bin-class}
\alias{ist.mod.bin}
\alias{ist.translator-class}
\alias{ist.translator}
\alias{ist.signatures-class}
\alias{ist.signatures}
\alias{ist.discriminator-class}
\alias{ist.discriminator}
\alias{ist.pathways-class}
\alias{ist.pathways}
\alias{ist.results-class}
\alias{ist.results}
\title{One-class model class}
\description{
One-class model class

Binary classifier model class

Translator class (fold changes and orthology mapping)

Signatures class (fold changes, metadata, orthology mapping)

Discriminator class (one-class/binary model and pathways)

Discriminator class (one-class/binary model and pathways)

Results class (one-class/binary pathways)
}
\section{Slots}{

\describe{
\item{\code{genes.mod}}{character vector, genes eligible to be in the model}

\item{\code{genes.top}}{character vector, genes in the model}

\item{\code{svm.tuned}}{\code{tune} object}

\item{\code{svm.best}}{\code{svm} object}

\item{\code{genes.mod}}{character vector, genes eligible to be in the model}

\item{\code{pls.mod}}{\code{mvr} object, from package \code{pls}}

\item{\code{Ncomp}}{numeric, optimal components for the PLS model}

\item{\code{train}}{\code{train} object, from package \code{caret}}

\item{\code{animal.data.table}}{data.frame with the animal fold changes}

\item{\code{org.from,org.to}}{character, organism to translate from and to,
should follow the ENSEMBL notation (hsapiens, mmusculus, etc)}

\item{\code{mapping.data.table}}{data.frame with the orthology mapping
between both organisms}

\item{\code{mapping.data}}{numeric named vector with the fold changes to
be added}

\item{\code{list.sig}}{list with data.frames containing the animal fold changes}

\item{\code{tab.meta}}{data.frame with the metadata - must contain the columns
\code{c("sig.id", "sig.name", "sig.org")}.}

\item{\code{org.to}}{character, organism to translate to,
should follow the ENSEMBL notation (hsapiens, mmusculus, etc)}

\item{\code{list.mapping}}{list with data.frames, by organism
(to translate from) name}

\item{\code{list.data}}{list with numeric named vectors with the fold changes to
be added}

\item{\code{X}}{numeric matrix, independent variables}

\item{\code{y}}{numeric vector, dependent variable}

\item{\code{org.to}}{character, organism with ENSEMBL notation
(hsapiens, mmusculus, etc)}

\item{\code{de.genes}}{character, disease-associated genes
(must be a subset of X's colnames)}

\item{\code{id.oc,id.bin,id.ref,id.ist}}{logical vectors, with length equal
to \code{nrow(X)}, indicating which samples should be used to fit
the one class, binary models, and be treated as reference and target
(in silico treatment should ideally move the reference to the target group)}

\item{\code{pathways.data.table}}{data.frame with pathway annotations}

\item{\code{iss,isd,isa}}{fitted models}

\item{\code{X}}{numeric matrix, independent variables}

\item{\code{y}}{numeric vector, dependent variable}

\item{\code{org.to}}{character, organism with ENSEMBL notation
(hsapiens, mmusculus, etc)}

\item{\code{pathways.table.oc,pathways.table.bin}}{data.frame with
pathway annotations for pathways expected to be unchanged/changed}

\item{\code{list.oc,list.bin}}{list with fitted models for
unchanged/changed pathways}

\item{\code{ist.signatures,ist.pathways}}{objects to generate the results from}

\item{\code{vec.gene2label}}{nicer labels for plotting genes}

\item{\code{group}}{factor, groups to stratify samples in plotting}

\item{\code{tab.decisions.oc,tab.metrics.oc}}{data.table objects for OC models}

\item{\code{tab.decisions.bin,tab.metrics.bin,tab.pls.bin}}{data.table objects
for BIN models}
}}

