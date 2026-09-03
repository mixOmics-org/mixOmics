
<!-- PLEASE DO NOT EDIT ./README.md BY HAND, EDIT ./inst/README.Rmd AND RENDER TO CREATE ./README.md -->

[![download](http://www.bioconductor.org/shields/downloads/release/mixOmics.svg)](https://bioconductor.org/packages/stats/bioc/mixOmics)
[![R build
status](https://github.com/mixOmics-org/mixOmics/workflows/R-CMD-check.yml/badge.svg)](https://github.com/mixOmics-org/mixOmics/actions)
[![](https://img.shields.io/github/last-commit/mixOmics-org/mixOmics.svg)](https://github.com/mixOmics-org/mixOmics/commits/master)
[![](https://codecov.io/gh/mixOmics-org/mixOmics/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mixOmics-org/mixOmics)
[![license](https://img.shields.io/badge/license-GPL%20(%3E=%202)-lightgrey.svg)](https://choosealicense.com/)
[![dependencies](http://bioconductor.org/shields/dependencies/release/mixOmics.svg)](http://bioconductor.org/packages/release/bioc/html/mixOmics.html#since)

![](https://raw.githubusercontent.com/mixOmics-org/.github/main/profile/assets/mixomics-logo.svg)

This repository contains the `R` package which is [hosted on
Bioconductor](http://bioconductor.org/packages/release/bioc/html/mixOmics.html)
and our development `GitHub` versions. Go to <https://mixomics.org> for
information on how to use mixOmics.

## Installation

(**macOS users only:** Ensure you have installed
[XQuartz](https://www.xquartz.org/) first.)

### From Bioconductor (recommended)

The best way to install `mixOmics` is using `Bioconductor`. You can see
the landing page for the release version of `mixOmics` on Bioconductor
[here](https://bioconductor.org/packages/release/bioc/html/mixOmics.html).
Make sure you have the latest R version and the latest `BiocManager`
package installed following [these
instructions](https://www.bioconductor.org/install/).

``` r
## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## install mixOmics
BiocManager::install('mixOmics')

## load mixOmics
library(mixOmics) 
```

### From Github

Bioconductor versions are updated twice a year, between these updates
you can download the latest version of `mixOmics` from `Github`. Note
that this latest version of mixOmics is under development and may not be
stable, check the gitHub page for releases which have passed package
testing.

``` r
## install devtools
install.packages("devtools")

## install latest github version of mixOmics
devtools::install_github("mixOmics-org/mixOmics")
```

### From Docker container

You can install our latest stable Github version of `mixOmics` via our
Docker container. You can do this by downloading and using the Docker
desktop application or via the command line as described below.

<details>

<summary>

Click to expand
</summary>

**Note: this requires root privileges**

1)  Install Docker following instructions at
    <https://docs.docker.com/docker-for-mac/install/>

**if your OS is not compatible with the latest version** download an
older version of Docker from the following link:

- MacOS: <https://docs.docker.com/docker-for-mac/release-notes/>
- Windows: <https://docs.docker.com/docker-for-windows/release-notes/>

Then open your system’s command line interface (e.g. Terminal for MacOS
and Command Promot for Windows) for the following steps.

**MacOS users only:** you will need to launch Docker Desktop to activate
your root privileges before running any docker commands from the command
line.

2)  Pull mixOmics container

``` bash
docker pull mixomicsteam/mixomics
```

3)  Ensure it is installed

The following command lists the running images:

``` bash
docker images
```

This lists the installed images. The output should be something similar
to the following:

    $ docker images 
      > REPOSITORY                       TAG       IMAGE ID       CREATED         SIZE
      > mixomicsteam/mixomics            latest    e755393ac247   2 weeks ago     4.38GB

4)  Activate the container

Running the following command activates the container. You must change
`your_password` to a custom password of your own. You can also customise
ports (8787:8787) if desired/necessary. see
<https://docs.docker.com/config/containers/container-networking/> for
details.

``` bash
docker run -e PASSWORD=your_password --rm -p 8787:8787 mixomicsteam/mixomics
```

5)  Run

In your web browser, go to `http://localhost:8787/` (change port if
necessary) and login with the following credentials:

*username*: rstudio\
*password*: (your_password set in step 4)

6)  Inspect/stop

The following command lists the running containers:

``` bash
sudo docker ps
```

The output should be something similar to the following:

``` bash
$ sudo docker ps
  > CONTAINER ID   IMAGE                   COMMAND   CREATED         STATUS         PORTS                    NAMES
  > f14b0bc28326   mixomicsteam/mixomics   "/init"   7 minutes ago   Up 7 minutes   0.0.0.0:8787->8787/tcp   compassionate_mestorf
```

The listed image ID can then be used to stop the container (here
`f14b0bc28326`)

``` bash
docker stop f14b0bc28326
```

</details>

## Contribution

We welcome community contributions concordant with [our code of
conduct](https://github.com/mixOmics-org/mixOmics/blob/master/CODE_OF_CONDUCT.md).
We strongly recommend adhering to [Bioconductor’s coding
guide](https://bioconductor.org/developers/how-to/coding-style/) for
software consistency if you wish to contribute to `mixOmics` R codes.

### Bug reports and pull requests

To report a bug (or offer a solution for a bug!) visit:
<https://github.com/mixOmicsTeam/mixOmics/issues>. We fully welcome and
appreciate well-formatted and detailed pull requests. Preferably with
tests on our datasets.

<details>

<summary>

Set up development environment
</summary>

- Install the latest version of R
- Install RStudio
- Clone this repo, checkout master branch, pull origin and then run:

``` r
install.packages("renv", Ncpus=4)
install.packages("devtools", Ncpus=4)

# restore the renv environment
renv::restore()

# or to initialise renv
# renv::init(bioconductor = TRUE)

# update the renv environment if needed
# renv::snapshot()

# test installation
devtools::install()
devtools::test()

# complete package check (takes a while)
devtools::check()
```

</details>

### Discussion forum

We wish to make our discussions transparent so please direct your
analysis questions to our discussion forum
<https://mixomics-users.discourse.group>. This forum is aimed to host
discussions on choices of multivariate analyses, as well as comments and
suggestions to improve the package. We hope to create an active
community of users, data analysts, developers and R programmers alike!
Thank you!

## What’s New

#### May 2026

- bug fix implemented for
  [\#377](https://github.com/mixOmicsTeam/mixOmics/issues/377) fixed
  vector `xlim` handling in the `graphics` branch of `plotLoadings()`
  contribution plots
- enhancement request
  [\#381](https://github.com/mixOmicsTeam/mixOmics/issues/381) expanded
  test coverage for model utilities and internal helper paths, with
  small fixes made for issues exposed by the new tests

#### April 2026

\*\* Version 6.36.0 \*\*

[Bioconductor release version
6.36.0](https://www.bioconductor.org/packages/release/bioc/html/mixOmics.html)
released 29th April 2026 mixOmics can now be downloaded using
[Bioconductor version
3.23](https://bioconductor.org/news/bioc_3_23_release/) and is
compatible with R 4.6.0

- bug fix implemented for
  [\#374](https://github.com/mixOmicsTeam/mixOmics/issues/374) replaced
  deprecated `aes_string()` with `aes(.data[[...]])` in plotting
  functions
- bug fix implemented for
  [\#379](https://github.com/mixOmicsTeam/mixOmics/issues/379) replaced
  deprecated `size` aesthetic with `linewidth` in ggplot2 line geoms

#### April 2025

\*\* Version 6.32.0 \*\*

[Bioconductor release version
6.32.0](https://www.bioconductor.org/packages/release/bioc/html/mixOmics.html)
released 17th April 2025 mixOmics can now be downloaded using
[Bioconductor version
3.21](https://bioconductor.org/news/bioc_3_21_release/) and is
compatible with R 4.5.0

- feature request
  [\#345](https://github.com/mixOmicsTeam/mixOmics/issues/345) updated
  functionality for `plotLoadings()` so can plot in ggplot2 style and
  customise aesthetics
- bug fix implemented for
  [\#357](https://github.com/mixOmicsTeam/mixOmics/issues/357)
  `plotIndiv()` not handing `pch` ordering correctly
- enhancement request
  [\#332](https://github.com/mixOmicsTeam/mixOmics/issues/332) increased
  test coverage for `plotIndiv()`
- implemented a new unit testing framework for plotting functions using
  `vdiffr` package

#### March 2025

- enhancement request
  [\#353](https://github.com/mixOmicsTeam/mixOmics/issues/353) better
  error message in `perf()` when one sample in a class
- enhancement request
  [\#340](https://github.com/mixOmicsTeam/mixOmics/issues/340) expand
  test coverage for main functions
- enhancement request
  [\#336](https://github.com/mixOmicsTeam/mixOmics/issues/336)
  streamline multiblock functions by removing `scheme` and `init` args

Also explored potential unusual behaviour of: \* zero variance handling
in `block.splsda`
[\#352](https://github.com/mixOmicsTeam/mixOmics/issues/352) \* `perf()`
giving non-intuitve per-class error rates
[\#355](https://github.com/mixOmicsTeam/mixOmics/issues/355)

#### November 2024

- enhancement request
  [\#216](https://github.com/mixOmicsTeam/mixOmics/issues/216)
  implemented parallel processing using `BPPARAM` across all `tune()`
  functions
- feature request
  [\#335](https://github.com/mixOmicsTeam/mixOmics/issues/335) added
  `seed` argument to `perf()` functions for better reproducibility
- feature request
  [\#334](https://github.com/mixOmicsTeam/mixOmics/issues/334) added
  `seed` argument to `tune()` functions for better reproducibility
- bug fix implemented for
  [\#303](https://github.com/mixOmicsTeam/mixOmics/issues/303) multiple
  solutions found in `perf()` returns error
- bug fix implemented for
  [\#307](https://github.com/mixOmicsTeam/mixOmics/issues/307)
  `plotIndiv()` ellipses colours not matching points, now sample group
  order is respected and colours can be customised for points and
  ellipses
- updated documentation to fix issue
  [\#297](https://github.com/mixOmicsTeam/mixOmics/issues/297) broken
  link in bookdown vignette
- updated documentation to fix issue
  [\#296](https://github.com/mixOmicsTeam/mixOmics/issues/296) typo in
  vignette

The performance assessment and parameter tuning workflow has been
streamlined as described in issue
[\#343](https://github.com/mixOmicsTeam/mixOmics/issues/343)

- New function: `perf.assess()` This function essentially runs `perf()`
  on final model but only returns performance metrics for the number of
  components used in the final model. Designed to be used in the final
  step of the workflow for quantifying final model performance. Outputs
  a list of values but no plotting functionality avaliable. See [PR
  \#344](https://github.com/mixOmicsTeam/mixOmics/pull/344) for more
  details.

- Additional functionality for `tune()` functions and new `tune()`
  functions created `tune()` can now be used in its original capacity
  (to tune number of variables and components simultaneously) or just to
  tune number of components by internally calling `perf()`. Designed to
  be used for tuning both components and variables to keep across
  (s)PCA, (s)PLS, (s)PLSDA, block (s)PLSDA and mint (s)PLSDA models See
  [PR \#348](https://github.com/mixOmicsTeam/mixOmics/pull/348) for more
  details.

#### October 2024

\*\* Version 6.30.0 \*\*

[Bioconductor release version
6.30.0](https://bioconductor.org/packages/release/bioc/html/mixOmics.html)
released end of October 2024 Minor bug fixes and updated deprecated code
and unit tests, no major code changes and no changes to user experience
of mixOmics.

- bug fix implemented for
  [\#293](https://github.com/mixOmicsTeam/mixOmics/issues/293)
  `splsda()` example code error

#### March 2022

- bug fix implemented for [Issue
  \#196](https://github.com/mixOmicsTeam/mixOmics/issues/196). `perf()`
  can now handle features with a `(s)pls` which have near zero variance.
- bug fix implemented for [Issue
  \#192](https://github.com/mixOmicsTeam/mixOmics/issues/192).
  `predict()` can now handle when the testing and training data have
  their columns in different orders.
- bug fix implemented for [Issue
  \#178](https://github.com/mixOmicsTeam/mixOmics/issues/178). If the
  `indY` parameter is used in `block.spls()`, `circosPlot()` can now
  properly identify the $Y$ dataframe.
- bug fix implemented for [Issue
  \#172](https://github.com/mixOmicsTeam/mixOmics/issues/172). `perf()`
  now returns values for the `choice.ncomp` component when `nrepeat`
  $< 3$ whereas before it would just return `NA`s.
- bug fix implemented for [Issue
  \#171](https://github.com/mixOmicsTeam/mixOmics/issues/171). `cim()`
  now can take `pca` objects as input.
- bug fix implemented for [Issue
  \#161](https://github.com/mixOmicsTeam/mixOmics/issues/161).
  `tune.spca()` can now handle `NA` values appropriately.
- bug fix implemented for [Issue
  \#150](https://github.com/mixOmicsTeam/mixOmics/issues/150). Provided
  users with a specific error message for when `plotArrow()` is run on a
  `(mint).(s)plsda` object.
- bug fix implemented for [Issue
  \#122](https://github.com/mixOmicsTeam/mixOmics/issues/122). Provided
  users with a specific error message for when a `splsda` object that
  has only one sample associated with a given class is passed to
  `perf()`.
- bug fix implemented for [Issue
  \#120](https://github.com/mixOmicsTeam/mixOmics/issues/120).
  `plotLoadings()` now returns the loading values for features from
  **all** dataframes rather than just the last one when operating on a
  `(mint).(block).(s)plsda` object.
- bug fix implemented for [Issue
  \#43](https://github.com/mixOmicsTeam/mixOmics/issues/43). Homogenised
  the way in which `tune.mint.splsda()` and `perf.mint.splsda()`
  calculate balanced error rate (BER) as there was disparity between
  them. Also made the global BER a weighted average of BERs across each
  study.
- enhancement implemented for [Issue
  \#30/#34](https://github.com/mixOmicsTeam/mixOmics/issues/34). The
  parameter `verbose.call` was added to most of the methods. This
  parameter allows users to access the specific values input into the
  call of a function from its output.
- bug fix implemented for [Issue
  \#24](https://github.com/mixOmicsTeam/mixOmics/issues/24).
  `background.predict()` can now operate on `mint.splsda` objects and
  can be used as part of `plotIndiv()`.

#### July 2021

- new function `plotMarkers` to visualise the selected features in block
  analyses (see <https://github.com/mixOmicsTeam/mixOmics/issues/134>)
- `tune.spls` now able to tune the selected variables on both `X` and
  `Y`. See `?tune.spls`
- new function `impute.nipals` to impute missing values using the nipals
  algorithm
- new function `tune.spca` to tune the number of selected variables for
  pca components
- `circosPlot` now has methods for `block.spls` objects. It can now
  handle similar feature names across blocks. It is also much more
  customisable. See advanced arguments in `?circosPlot`
- new `biplot` function for `pca` and `pls` objects. See
  `?mixOmics::biplot`
- `plotDiablo` now takes `col.per.group` (see \#119)

#### April 2020

- weighted consensus plots for DIABLO objects now consider per-component
  weights

#### March 2020

- `plotIndiv` now supports (weighted) consensus plots for block
  analyses. See the example in [this
  issue](https://github.com/mixOmicsTeam/mixOmics/issues/57)
- `plotIndiv(..., ind.names=FALSE)` [warning
  issue](https://github.com/mixOmicsTeam/mixOmics/issues/59) now fixed

#### January 2020

- `perf.block.splsda` now supports calculation of combined AUC
- `block.splsda` bug which could drop some classes with
  `near.zero.variance=TRUE` now fixed
