% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit-ist}
\alias{fit-ist}
\alias{fit.ist.translator}
\alias{fit.ist.signatures}
\alias{fit.ist.discriminator}
\alias{fit.ist.pathways}
\alias{fit.ist.oc}
\alias{fit.ist.bin}
\alias{fit.ist.results}
\title{Fit methods}
\source{
https://stackoverflow.com/questions/35907477/caret-package-
stratified-cross-validation-in-train-function
}
\usage{
fit.ist.translator(object, ...)

fit.ist.signatures(object, ...)

fit.ist.discriminator(object, flavour, ist.translator, ...)

fit.ist.pathways(object, ...)

fit.ist.oc(
  X.ref,
  de.genes,
  n.genes = 500,
  seed = 1,
  args.oc.svm = default.args.oc.svm(),
  args.oc.train = default.args.oc.train(),
  ...
)

fit.ist.bin(
  X,
  y,
  de.genes,
  seed = 1,
  args.bin.train = default.args.bin.train(),
  args.bin.pls = default.args.bin.pls(),
  ...
)

fit.ist.results(object, ...)
}
\arguments{
\item{object}{an \code{ist.translator} or an
\code{ist.discriminator} object}

\item{...}{in \code{fit.ist.discriminator()},
these parameters are passed to
\code{fit.mod.oc()}, \code{fit.mod.bin()} and optionally
to \code{get.de.genes()}.
Ignored otherwise}

\item{flavour}{flavour of in silico treatment;
can contain iss (one-class) and/or
isd, isa (binary classifier)}

\item{ist.translator}{an \code{ist.translator} object (for isa),
to obtain the list of genes allowed in the model}

\item{X.ref}{a matrix containing the target organism and group for iss}

\item{de.genes}{disease genes, which will be EXCLUDED from the
training in iss, USED in isd and IGNORED in isa}

\item{n.genes}{numeric, maximum number of genes to include in the iss model.
Currently genes with the maximum mean expression are kept}

\item{seed}{integer, seed before train and final model fit, or \code{NA}
not to set any}

\item{X}{matrix with the regressors for isd/isa}

\item{y}{numeric response for isd/isa}

\item{args.bin.train}{list with options for
\code{caret::train} in isd/isa}

\item{args.bin.pls}{list with options for \code{pls::plsr} in isd/isa}

\item{gamma, nu, kernel, tunecontrol}{parameters for iss, passed to
\code{e1071::tune.svm}}
}
\value{
a fitted object with class \code{ist.translator},
\code{ist.discriminator}, \code{ist.mod.oc} or \code{ist.mod.bin}
}
\description{
\code{fit.ist.translator()} fits the translator by retrieving
orthologs and intersecting the identifiers with those in the fold changes.
If the translator was initialised with an ortholog mapping,
it will be reused, this can avoid the slow, repeated querying of the same
orthology mapping when dealing with multiple translators.

\code{fit.ist.signatures()} fits the translator on multiple
signatures by first stratifying by organism, retrieving the
orthologs and intersecting the identifiers with those in the fold changes.
A custom orthology mapping can be used.

\code{fit.ist.discriminator()} fits the
discriminator on the org.to
organism

\code{fit.ist.pathways()} fits the
pathway models on the org.to organism (must have pathway data),
by first intersecting the pathways and then fitting the oc/bin models

\code{fit.ist.oc()} fits the one-class
support vector machine

\code{fit.ist.bin()} fits the binary classifier,
a PLS regression

\code{fit.ist.results()} applies a battery of signatures
in an \code{ist.signatures} object to the pathway models in
\code{ist.pathways}. Given the complexity of the results, an object of
class \code{ist.results} is returned.
}
