#' Omics Data Integration Project
#'
#' Multivariate methods are well suited to large omics data sets where the
#' number of variables (e.g. genes, proteins, metabolites) is much larger than
#' the number of samples (patients, cells, mice). They have the appealing
#' properties of reducing the dimension of the data by using instrumental
#' variables (components), which are defined as combinations of all variables.
#' Those components are then used to produce useful graphical outputs that
#' enable better understanding of the relationships and correlation structures
#' between the different data sets that are integrated.
#'
#' mixOmics offers a wide
#' range of multivariate methods for the exploration and integration of
#' biological datasets with a particular focus on variable selection. The
#' package proposes several sparse multivariate models we have developed to
#' identify the key variables that are highly correlated, and/or explain the
#' biological outcome of interest. The data that can be analysed with mixOmics
#' may come from high throughput sequencing technologies, such as omics data
#' (transcriptomics, metabolomics, proteomics, metagenomics etc) but also beyond
#' the realm of omics (e.g. spectral imaging).
#'
#' The methods implemented in
#' mixOmics can also handle missing values without having to delete entire rows
#' with missing data. A non exhaustive list of methods include variants of
#' generalised Canonical Correlation Analysis, sparse Partial Least Squares and
#' sparse Discriminant Analysis. Recently we implemented integrative methods to
#' combine multiple data sets: N-integration with variants of Generalised
#' Canonical Correlation Analysis and P-integration with variants of multi-group
#' Partial Least Squares.
#'
#' @aliases mixOmics-package
#' @import MASS lattice igraph ggplot2 corpcor parallel RColorBrewer
#' @importFrom grDevices as.graphicsAnnot chull col2rgb colorRamp colorRampPalette colors dev.cur dev.new dev.off dev.prev dev.set devAskNewPage graphics.off gray gray.colors heat.colors rgb jpeg pdf tiff x11 adjustcolor rainbow
#' @importFrom graphics abline arrows axis barplot box image layout legend lines locator mtext par plot plot.default points polygon rect segments strheight strwidth symbols text title Axis boxplot rasterImage matplot
#' @importFrom stats cor cov dist hclust lm lsfit median na.omit order.dendrogram predict quantile reorder var sd pnorm aggregate t.test
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom dplyr arrange rename filter group_by mutate n row_number summarise ungroup
#' @importFrom rlang .data
"_PACKAGE"
