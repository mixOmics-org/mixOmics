## version 6.38.0

Bioconductor 3.24.

### new features / enhancements / changes

* mixOmics is now licensed under AGPL (>= 3) instead of GPL (>= 2) (#421)
* `circosPlot` is drawn by new code (#419)
* `predict` on multi-block models is 3 to 67 times faster, and `perf` and `tune.block.plsda` about 2 times faster, with identical results (#369, #383)
* `perf` and `perf.assess` on DIABLO models accept `folds` as a list (#400)
* `mint.block.pls`, `mint.block.spls`, `mint.block.plsda` and `mint.block.splsda` return model objects with the same structure as the `block.*` models, so `predict` works on them (#386)
* `?plot.perf.pls` now has examples (#414)
* `citation("mixOmics")` gives the full author names and the DOI (#416)
* the GitHub repository moved from `mixOmicsTeam/mixOmics` to `mixOmics-org/mixOmics`
* expanded test coverage for model utilities and internal helpers (#381)

### bug fixes

* `circosPlot` no longer fails when a block, including the Y block, has a single selected variable (#402, #324)
* `circosPlot` showed block names containing "chr" truncated and failed with `line = TRUE` (#419)
* `circosPlot` with several blocks in `blocks.link` kept only about half of the links (#419)
* `plotVar` no longer places points outside the correlation circle when the plotted components are not orthogonal (#405)
* `keepX` and `keepY` of `block.spls` and `mint.block.*` models now give the values used (#386)
* `plotLoadings` handles a vector `xlim` for contribution plots in the `graphics` style (#377)

## version 6.36.0

Bioconductor 3.23, April 2026.

### new features / enhancements / changes

* plotting functions use `aes(.data[[...]])` instead of `aes_string()`, which is deprecated in ggplot2 (#374)
* ggplot2 line geoms use `linewidth` instead of the deprecated `size` (#379)
* gsignal is no longer a dependency

### bug fixes

* `plotIndiv` on MINT models now plots a single study selected with `study` (#364)

## version 6.34.0

Bioconductor 3.22, October 2025. No changes from 6.32.0.

## version 6.32.0

Bioconductor 3.21, April 2025.

### new features / enhancements / changes

* new function `perf.assess` to assess the performance of a final model, for the number of components it uses (#343)
* `tune` can tune the number of components only, by calling `perf`, or the number of components and variables together, for (s)PCA, (s)PLS, (s)PLS-DA, block (s)PLS-DA and MINT (s)PLS-DA (#343)
* all `tune` functions support parallel processing through `BPPARAM` (#216)
* `perf` and `tune` functions take a `seed` argument for reproducibility (#334, #335)
* `plotLoadings` can plot in ggplot2 style with customisable aesthetics (#345)
* the `scheme` and `init` arguments of the multi-block functions are removed (#336)
* `perf` gives a clearer error when a class has a single sample (#353)
* expanded test coverage, including visual tests of the plots with vdiffr (#332, #340)

### bug fixes

* `plotIndiv` ellipses now have the same colours as the points, follow the order of the groups, and can be coloured separately (#307)
* `plotIndiv` uses `pch` in the right order (#357)
* broken link to the bookdown vignette (#297) and typo in `?summary.mixo_pls` (#296) fixed

## version 6.30.0

Bioconductor 3.20, October 2024. Dependencies updated and BiocCheck errors fixed; no change for users.

### bug fixes

* `splsda` and `plotLoadings` help page examples fixed (#293)

## version 6.28.0

Bioconductor 3.19, April 2024. No changes from 6.26.0.

## version 6.26.0

Bioconductor 3.18, October 2023. Fixed the vignette build on Bioconductor; no change to the package code.

## version 6.24.0

Bioconductor 3.17, April 2023.

### new features / enhancements / changes

* `auroc` takes a list of (s)PLS-DA models and draws their ROC curves in one plot (#242)
* `auroc` works on `mint.block.(s)plsda` models (#246)
* `circosPlot` works on `block.pls` and `block.plsda` models (#248)
* `cim` works on `block.spls` models (#187)
* `pca` output now includes `x` and `rotation`, as `prcomp` does, so `cim` accepts `pca` objects (#171)
* `network` takes `size.node` and `graph.scale` to set the size of the nodes and of the plot (#233)
* `plotIndiv` on MINT models can show the prediction areas from `background.predict` (#24)
* `tune.spls` supports parallel processing (#214)
* main functions take `verbose.call` to return the call with all argument values (#30, #34)
* clearer errors when Y has missing values in DA methods (#111), when a class has a single sample in `perf` (#122), and in `tune.spca` (#275)

### bug fixes

* `perf.mint.splsda` and `tune.mint.splsda` give consistent error rates (#43)
* `perf` on `plsda` no longer returns `NA` for `choice.ncomp` (#172)
* `perf.pls` error "non-conformable matrices" fixed (#286)
* `plotArrow` error "Aesthetics must be either length 1 or the same as the data" fixed (#150)
* online vignette updated (#230)

## version 6.22.0

Bioconductor 3.16, November 2022.

### new features / enhancements / changes

* `network` takes `plot.graph` to return the network without plotting it (#145)
* `circosPlot` takes `legend.title` (#260)
* `plotLoadings` returns the loadings of all blocks (#120)

### bug fixes

* `circosPlot` works on `block.spls` models fitted with `indY` (#178)
* `predict` on DIABLO models handles test data whose variables are in a different order (#192)
* `predict` no longer returns `NaN` for a whole block (#266) and handles near-zero-variance features (#268)
* `perf` on `pls` handles near-zero-variance features (#196)
* `tune.block.splsda` with leave-one-out validation fixed (#213)
* `plotVar` works on multi-block (s)PLS1 models (#237)
* documentation of `perf.pls` (#130) and `plot.tune` (#189) fixed

## version 6.20.0

Bioconductor 3.15, April 2022.

### bug fixes

* `tune.spls` selected the wrong `keepX` and `keepY` because of an incorrect t-test alternative hypothesis (#164)
* `plot` on `tune.mint.splsda` objects fixed, and the `sd` argument of `plot.tune` restored (#167)
* compatibility with R 4.2 (#162, #208)

## version 6.18.0

### new features / enhancements / changes

* new function `plotMarkers` to visualise the selected features in block analyses (see #134)
* `auroc` title now fixed (#135)
* `cimDiablo` takes `trim` argument to customise outlier filtering (#136)
* `plotIndiv.pca` default shape set to `16`
* `circosPlot` & `network` now support blocks with similar feature names 
* `circosPlot` now has methods for block.spls objects
* `circosPlot` now takes new formal and advanced args for customisation. See `?circosPlot`.

### bug fixes

* `plotVar` legend colour mismatch bug fixed
* `plotDiablo` error `undefined variable (Y)` fixed
* `nipals` initialisation bug with high-variance high-NA rate column fixed 
* `cim` bug with high NA rate data fixed using imputation by the column mean

## version 6.16.0

### new features / enhancements / changes

* `biplot` default colours now match `plotIndiv`
* `biplot` now takes `cex` argument
* `nipals` now takes `center` and `scale` arguments
* `nipals` now only outputs `p`, `t` and `eig`
* new function `impute.nipals` to impute missing values using NIPALS algorithm
* `nipals` function checks for orthogonality of components when NA present
* `plot.tune` legends improved
* `plot.tune` now uses colour-blind friendly colours
* new `tune.spls` function to perform variable selection on both X and Y
* `tune.spls` now chooses optimal keepX even if nrepeat < 3
* `tune.spls` now takes `validation` arg
* `tune.spca` is now much faster
* `circosPlot` links can now have adjustable width (#118)
* `plotDiablo` now takes `col.per.group` (#119)
* terminology change: `consensus` renamed to `average` in Diablo context

### bug fixes

* `plotVar` bug fixed

## version 6.14.0

### new features / enhancements / changes

* `circosPlot`: The radial location of feature names can now be cutomised using `var.adj`
* added `plot` and `print` methods for `nipals` ouput (#87)
* all Discriminat Analyses now run solely on `mode=regression` (#79)
* `cim` argument change: `threshold` replaced by `cutoff`
* `nipals` and `pca` with missing values allow skipping reconstitution of the input matrix
* `tune.block.splsda` now allows random number seed also for parallel processing (#72)
* New `biplot` methods for the pca family (#90)

### bug fixes

* `plotIndiv`: Legend bug which misspecified the groups resolved
* `plotIndiv`: Legends now ordered as inputted, and not alphabetically
* `plot` method issue for `spca` resolved 
* `plotLoadings.spca` bug with `var.names` now fixed (#81)
* `ipca` deprecation warning fixed

-------------------------------------------------------------------------------
## version 6.12.0

### new features / enhancements

* `plotLoadings`'s infamous *figure margins too large* error now handled and informative condition thrown
* `circosPlot`'s `lines` argument default to `FALSE` now
* `circosPlot`'s inconsistentcy of blocks with identical `X` names fixed
*  consensus and weighted consensus plots now supported for `plotIndiv` with relevant block analyses
* `plotLoadings`'s feature name trimming can be customised
* `block.splsda` bug which could drop some Y factors with `near.zero.variance=TRUE` fixed
* `perf.block.splsda` now supports calculation of combined and per-block AUC
* model improvement significance can be custmoised in all `perf` and `tune` functions
* `perf.block.splsda` is now much faster and supports FORK clusters
* `tune.(s)pls(da)`, `perf.(s)plsda` now support FORK clusters

### bug fixes

* `circosPlot`'s faded `lines` bug when many `NA`s present fixed
* `tune.block.splsda()` bug when using fixed `test.keepX` over two or more blocks fixed
* `circosPlot` and `plotLoadings` bug caused by features with `NA`s fixed
* `plotIndiv(..., ind.names = FALSE)` warning/bug fixed
* `tune.block.splsda` bug on Windows parallelisation fixed
* `perf` and `tune` functions' issue  when choosing the optimum component resolved
* added option to suppress `auroc` from printing all the AUCs
-------------------------------------------------------------------------------
## version 6.10.0

### new features / enhancements

* parallel processing on `tune.block.splsda` improved
* `tune.block.splsda` now supports more distances
* You can now customise `auroc` plots. Refer to documentation for more info

### bug fixes

* single factor multilevel error in `pls` fixed
* fixed over-estimated correlation of `cim` for `mixo_(s)pls` objects with single component 
* margin error in `cim` now handled properly
* fixed `plotLoadings` error for very long variable names
* `predict` function bug for single sample prediction fixed
* `plotLoadings` bug for long variable names fixed
* Fixed `tune.spls` and `pef.plsda` bugs when using `cpus` argument for parallel 
processing
* `perf.plot` bug in extracting names fixed
* Few fixes for `tune.splsda` with AUC

### minor improvements

* missing values in `plotIndiv`'s `group` argument no more throws error
* `mixOmics::predict` function documentation now more accessible
* names of `linnerud` datasets fixed
* `plot.perf` now respects `ylim` arguments for custom y range
* package startup message with direct liks to useful resources
* `mixOmics` function documentation disambiguated with instruction on how to get
package help
* Updated onLoad message with discussion forum info, bug reports, and more
* Dropped legacy `comp.tol` argument from `pca`
* `plot.perf` now respects `ylim` arguments for custom y range
* Added Code of Conduct

-------------------------------------------------------------------------------
## version 6.8.0

* NOW HOSTED ON BIOCONDUCTOR

-------------------------------------------------------------------------------

## version < 6.8.0

* Refer to `./inst/legacy/NEWS-old` on GitHub repo