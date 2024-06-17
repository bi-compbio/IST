# NEWS file

### IST 0.5.30

* Fixed `NAMESPACE` exports for data objects
* Added more readable default labels

### IST 0.5.29

* Fixed bug where the vignette could not be built

### IST 0.5.28

* Fixed bug where `max.genes` would not work for ggplot/long table format

### IST 0.5.27

* Added the `max.genes` option to plot genemaps -- by default it keeps all genes (so nothing changes from prior versions), but allows limiting their number

### IST 0.5.26

* Minor fixes in `show` methods
* Now pathway metadata can be specified (optional), and the automatic
columns `mod` (oc or bin) and `n.genes` (after mapping to human data)
are added
* Pathway metadata can be displayed in pathwaymaps by using `vars.meta.path`

### IST 0.5.25

* Updates in vignette, including graphical scheme
* Updated syntax 
* Added section on pathway/sig prioritisation

### IST 0.5.24

* Added prioritisers for pathways and signatures
* Fixed inexistent default argument in `get.tab.pathwaymap`
* Legacy batch functions are skipped in unit testing

### IST 0.5.23

* Important change: `group` now belongs to `ist.results` objects, not to 
`ist.pathways`. This avoids having to fit pathways again to change 
the groups.
* Error raises if `id.ref` or `id.ist` are all `FALSE`
* The orthology mapping and the ensembl-symbol mapping are now accessible 
with `::` (lazydata was not set in DESCRIPTION)
* `ist.path` can be defined with empty pathway tables for `oc` or `bin`
* Added getters for data, reponse, gene and pathway ids, sample groups 
and their levels

### IST 0.5.22

* Fixed bug in pheatmap with NULL title
* Updated NAMESPACE

### IST 0.5.21

* Allowed cartesian in data.table joins, as this would raise exceptions
in when going to 100k+ rows in delta tables
* When `sig.ids` has no overlap with a pathway, `get.tab.genemap()` would
raise an error - now returns NULL with a message
* Added option of plot titles, and genemaps have the pathway name by default
* Added helpers `save.pheatmap()`, `save.genemaps()` and `save.boxplots()` 
to provide a default way of exporting the basic plots

### IST 0.5.20

* Minor modifications so that CI goes through

### IST 0.5.19

* Fixed issue of missing gene labels in `ist.results`
* Added two label mapping vectors ready to use in `ist.results` instead of 
data.frame
* Expanded options for `pheatmap` through the `args.pheatmap` argument
and its defauls
* Added `getGroups()`, `getMetaSig()` and refactored some slot accessions
* By default, `group` is ordered now

### IST 0.5.18

* Fix order issue in boxplots
* Fixed default boxplot: now shows original data groups
* Fixed colour scale in heatmaps

### IST 0.5.17

* Updated vignette
* Exported plotting functions

### IST 0.5.16

* Small fixes to `do.isx()` (legacy) and to `plot.ist.boxplot()` (bug)

### IST 0.5.15

* Fixed error in default gene labels
* Now signatures can be subsetted in every plot type with `sig.ids`
* Removed assertions of strict subsets in signature/pathway ids

### IST 0.5.14

* Added default orthology mappings `data.list.orth` and gene labels 
for plotting `data.dt.genelabels`
* Added the slot `@vec.gene2label` to `ist.results`; if provided, 
gene labels are switched in the genemaps

### IST 0.5.13

* Fixed original signature metadata, that had extra rows
* Moved `id.ist` and `id.ref` to the `ist.results` object
* Updated vignette with current signature and pathway classes
* Updated CI dependencies

### IST 0.5.12

* Added more getters, esp `getSignatures()` and `getPathways()`
* Added plotting functions: `plot.ist.boxplots()`, 
`plot.ist.genemaps()`, `plot.ist.pathwaymaps()`
* The `group` field for reference data is optional. If omitted, 
grouping uses the different `y` values.

### IST 0.5.11

* Essential tables are now ready in `ist.results` objects
* Getter functions (e.g. `getTabSignatures()`) ease slot accession in 
`ist.results` objects
* Added helpers to generate the tables

### IST 0.5.10

* Added first getters to retrieve results tables (TODO: unit testing)
* Fixed minor bug in imports

### IST 0.5.9

* Added `ist.results` class for holding pathway results
* Added `fit()` and `show()` methods

### IST 0.5.8

* Prepared parallel implementations of fitting and predicting `ist.pathways` 
* However, **now the SerialParam backend is used** due to MulticoreParam 
hanging for no apparent reason (probably missing exports)

### IST 0.5.7

* `ist.pathways` should accept a `BiocParallel` backend now 
for `fit()` and `predict()`
* Added helper function with repeated logic for fitting `ist.pathways`

### IST 0.5.6

* Added a helper with the logic of applying fold changes to data matrices, 
as it is shared between some classes
* Added `predict()` to `ist.signatures`

### IST 0.5.5

* Added `fit()` and `show()` method for `ist.signatures`
* Now `orthIDcon()` is a unique interface for mapping within and between species
* Added `BiocParallel` dependency

### IST 0.5.4

* Small change in `sample.env.ist`
* Added new class `ist.signatures` to hold signatures from varios organisms

### IST 0.5.3

* Added new sample data `sample.env.ist`, more similar to what we are 
actually using IST on. Contains signatures from several organisms 
with metadata and pathways with metadata.

### IST 0.5.2

* Added defaults to OC models (as they were present for BIN)
* Now OC should accept other kernels as well (untested)
* Added a slot for models, storing the sample names used for training
* Decisions have new columns: `sample.intrain`, `sample.ref`, `sample.ist`

### IST 0.5.1

* Fixed strange bug in `NAMESPACE` (watch out for this in the future), 
where `do.isx` would be exported as a `S3` method called `do`.
Now roxygen generates the correct file again
* Minor changes

## IST 0.5.0

* Important API change: `ref.y` does not exist anymore. 
Now one can control the samples included in the models using `id.oc`, `id.bin`
and those that are IST'd with `id.ref` and `id.ist`
* TODO: make use of these in `plot()` and batch mode

### IST 0.4.15

* `reduce.reg()` can now wait for unfinished jobs and pass arguments to `loadRegistry()`

### IST 0.4.14

* The wrapper `do.isx()` does not return the models now, to save disk space and hopefully lower execution time in `reduceResults()`

### IST 0.4.13

* Updated the main `README.md` to give an overview of IST's purpose

### IST 0.4.12

* Fixed bug in `get.de.genes()` that would always restrict to 100
differential genes
* Documented `pathways.table.oc` and `pathways.table.bin` in wrapper

### IST 0.4.11

* Seems like the `0.4.10` version batch jobs completed successfully 
(isd, isa, iss)
* Fixed example for `define()`
* Batch mode seems to work for pathway models

Comments:

* If pathway models are going to be independent of the translator, 
they can be sped up by just fitting them once and then running all the 
translators on them. Might need another wrapper like `do.isx()`.

### IST 0.4.10

This version contains important changes on how the wrappers work.
I am trying to have a common interface to avoid duplicated code.
This is achieved by combining ellipsis and explicit arguments that
are passed on, and leveraging the common 
`define()`, `fit()` and `predict()` API.

* **It is unclear whether batch model is working yet**.
* Removed duplicated entries in docs
* Defined a way to create objects: `define()`. It is like `new()` but 
it can ignore mismatching arguments in ellipsis without throwing error.
* Now `do.isx()` can run both `ist.translator` and `ist.pathways`.

### IST 0.4.9

* Fixed unit testing (change in behaviour of `limma` helper)
* Fixed duplicated documented arguments
* Fixed uses of `1:n`
* Fixed tabs (multiple of 4)
* Removed lines of length >80 
* Updated documentation format
* Added `biocViews`
* Updated maintainer/author fields to be compliant
* `BiocCheck()` gives 1W (set.seed) and 5N

### IST 0.4.8

* Batch jobs can now run `isd`, `isa` ans `iss`
* Changed name of wrapper to run the three of them
* Added imports from `stats` and `utils`
* Helper to find DE genes will only return their ids by default,
`topTable` can be obtained with a flag
* Batch can run with and without list of `de.genes`
* Jackknifing is not run by default on pls models

### IST 0.4.7

* Added `statmod` package for robust `limma` estimates

### IST 0.4.6

* Added `limma` for CI/CD
* Now `fit()` triggers the DE genes analysis
* Added unit testing in workflow
* Fixed unit testing

### IST 0.4.5

* Added function to trigger `limma` analysis to find `de.genes`
* Moved unit testing of defaults to their own file

### IST 0.4.4

* Fix bug where 0-confidence orthologs would make it through
* Intra-species mapping now has same format as inter-species
* Cleaned up log files

### IST 0.4.3

* Small fix to skip batch tests

### IST 0.4.2

* Added the required packages `SummarizedExperiment`, `batchtools` to 
CI/CD
* Build ignores `reg-` and `Meta` directories

### IST 0.4.1

This commit contains the very basics for batch analysis, focused on 
one flavour for now (ISD).
This affects a good amount of files.

* Main function: wrapper `do.batch.isd()`
* Default SGE configuration files `default.config.batch()`
* Helper for reducing: `reduce.reg()` (UNTESTED) 
* Added batch configuration files for regular job and e-mail notification
(just works for my e-mail right now) in `inst/templates`
* Basic examples and nit testing
and basic submitter wrapping `do.isd()`
* `.gitignore` will ignore directories starting by `reg-`

Other changes

* Helpers now have a separate file from defaults

To do

* Check if batch version is actually using predefined orthologs
* Make sure no small files are left on the working directory (logs...)
* Test the reduce function

## IST 0.4.0

* Added `do.isd()` wrapper, documented and with example

### IST 0.3.14

* `fit.ist.translator()` should re-use orthology mapping if provided
* Added test case in workflow to cover this

### IST 0.3.13

* Removed `mapOrthologs()` function, was apparently a duplicate of 
`orthIDcon()` and was never called 
* Changed mistaken `utils::new()` to `methods::new()`

### IST 0.3.12

* Added the option to specify a `proxy` in the ortholog mapper `orthIDcon()`, 
checked the function arguments and added ellipsis

### IST 0.3.11

* Fixed `.gitlab-ci.yml`

### IST 0.3.10

* Still dealing with vignette prebuilt index
* Info [in github](https://github.com/r-lib/devtools/issues/587) and 
in SO ([link1](https://stackoverflow.com/questions/39649133/package-has-a-vignettebuilder-field-but-no-prebuilt-vignette-index) [link2](https://stackoverflow.com/questions/30976308/how-do-i-prebuild-a-vignette-index-for-an-r-package))
* Now I think this warning is because of the `--no-build-vignettes` option
* Modified the `.gitlab-ci.yml` file because the error 
said `BiocStyle` was not available
* Importing `utils` because now new objects are created with `utils::new()`

### IST 0.3.9

* Trying to fix the following error in builder:
`Package has a VignetteBuilder field but no prebuilt vignette index.`
* Changed `.Rbuildignore` to not ignore `^doc$` and `^Meta$`

### IST 0.3.8

* Added quickstart vignette by taking the example from `data-raw/01_showcase.R`

### IST 0.3.7

* Fixed small erorr in showcase
* Small changes in `.Rbuildignore`, `.gitignore`
* Now only 20 repeats in the binary model training (down from 100)

### IST 0.3.6

* Expanded example scripts to try batch jobs

### IST 0.3.5

* Expanded example data for testing purposes

### IST 0.3.4

* Fixed small bug for intra-species translators

### IST 0.3.3

* Added `show` method for `ist.pathways`

### IST 0.3.2

* Added `predict` method for `ist.pathways`

### IST 0.3.1

* Added `fit` for `ist.pathways`
* Added pathways to the sample data

## IST 0.3.0

* Added `ist.pathways` class

### IST 0.2.4

* Small fix

### IST 0.2.3

* CI enabled

### IST 0.2.2

* The one-class model sets seeds as well
* Seeds are optional now (can be useful for parallel jobs)

### IST 0.2.1

* Added helper functions to generate default parameters for binary classifiers, 
was particularly annoying to set up seeds in `caret::train()`
* Removed internal data with defaults
* Changed `preProc` to `preProcess`, it's unclear why this is referred to as 
both in the documentation (`?train`)
* Seeds are controlled in the binary classifier, two approaches:
with a `seed` argument, which calls to `set.seed()` before training, 
and using the `seeds` built-in argument in `trainControl`. 
Still have to check how this works in a parallel setup.

## IST 0.2.0

* Changed the validation scheme for the binary classifier
* `plsr` - now uses `caret` and runs by default 100 rounds of repeated CV
* Extra slot with the `train` object
* Supress only the warning on regression using two values only
* When fitting the final PLS, all the specified number of components 
will be fit (especially for scoreplot purposes); however, when predicting 
the decision, only the optimal `Ncomp` will be taken into account
* By default, jackknifing is enabled on the final model

### IST 0.1.23

* Fixed bug in `predict.ist.bin()` that would not predict scores even
if asked to do so
* Added requirement on the checkmate package version

### IST 0.1.22

* Added `plot` method for scores and decision values
* Doc and test updates

### IST 0.1.21

* Added the `ref.y` slot to iss to have traceability
* Improved `show`

### IST 0.1.20

* Added true fold changes to sample data
* Showcase is longer and has density plots
* Improved documentation on `fit` and `predict` methods
* Added `show` methods

### IST 0.1.19

* Now `predict` can give the scores of the `pls` models (isd/isa)
* Defaults to 2 components
* Added example (future vignette)

### IST 0.1.18

* Refactored `type` to `flavour`

### IST 0.1.17

* Small class refactoring
* Fixed bug in `predict`

### IST 0.1.16

* Removed unnecessary file for `predict` method

### IST 0.1.15

* Added `predict` method for the discriminator

### IST 0.1.14

* Removed `isc`, now it is a particular case of `isa`
* Changes to the class definitions (more checks, slot removal)
* Primitive workflow already available
* `fit` method for the discriminator that applies to all flavours

### IST 0.1.13

* Refactoring of class and methods names (were too long)

### IST 0.1.12

* Added a `predict` method for binary classifier

### IST 0.1.11

* Added class for the classifier model in ISS
* Added a `fit` method for it

### IST 0.1.10

* Added a `predict` method for one-class

### IST 0.1.9

* Added class for the one-class model in ISS
* Added a `fit` method for it

### IST 0.1.8

* Same as version `0.1.7` with the predict method

### IST 0.1.7

* Now the fit method is in its own file and better commented

### IST 0.1.6

* Added predict method for the translator
* Modified the fit method, to join the mapping data.tables once only
* Updated the format of the sample data
* Updated unit testing

### IST 0.1.5

* Fixed fold change data frame (gene ids must be character)
* Added fit method for the translator

### IST 0.1.4

* Added intra-organism mapper

### IST 0.1.3

* Added classes for translator and discriminator

### IST 0.1.2

* Added sample data 
* Fixed doc

### IST 0.1.1

* Modified `DESCRIPTION`
* Added basic orthology mapper

## IST 0.1.0

* Initial addition of functions and tests
