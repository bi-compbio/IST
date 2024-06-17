# IST

IST stands for In Silico Treatment, a statistical approach for a data-driven selection of animal models and treatments. The main purpose of IST is to assess the actual translation of animal models and the potential of treatments with promising effects on murine models. The methodology is not limited to mouse, but can use any other organism with reasonable orthology mapping to humans, or even other human data itself.

## Input data

-   A human dataset with molecular readouts for control and disease (default: gene expression)
-   An animal model dataset, with the same types of molecular readouts; having the significant fold changes and their magnitude is enough
-   Gene sets / pathways to restrict translatability
-   A mapping between species (default: gene orthologs)

## Algorithm and classes

-   Fit `ist.signatures` object: define a way to translate changes between species. By default, assumes that the fold changes in humans translate directly to their ortholog animal genes.
-   Fit `ist.pathways` object: build statistical models that discriminate healthy from patients in human data, by gene set.
-   Fit `ist.results` object, using the two above: simulates animal model effect on human controls (or treatment effect in patients), then applies statistical models and quantifies direction and magnitude of change.

## Output

-   Decision values that quantify the faithfulness of the animal model signature to reproduce the molecular changes in human disease, or the effectiveness of a drug if it was applied to humans
-   Recapitulation (100% ideal, 0% no effect) of the human changes as observed in the animal data, per gene set / pathway
-   Contribution to the overall recapitulation by gene (%, positive indicates right direction whereas negative indicates opposite direction versus human data)

## Methods

Classes have the methods `fit()` and `predict()`, automatically dispatched according to the object class. See the `quickstart` vignette for a reproducible minimal example. Running `browseVignettes("IST")` will display all the available vignettes. Small synthetic exemplary data is bundled as well, see `?sample.data.ist`.

## Other repositories

This publication bundles three repositories:

* `IST`: contains the logic to apply In Silico Treatment
* `ISTResults`: use cases on Idiopathic Pulmonary Fibrosis (IPF) and Non-Alcoholic SteatoHepatitis (NASH)
* `ISTBrowser`: self-service shiny app to interactively browse the results of IST
