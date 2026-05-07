ARG BIOC_VERSION=devel
FROM bioconductor/bioconductor_docker:$BIOC_VERSION

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "for (p in c('remotes', 'devtools')) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, Ncpus = 2); Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE); options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::repositories(); remotes::install_deps(dependencies = TRUE, build_vignettes = FALSE, upgrade = 'never', repos = BiocManager::repositories()); devtools::install('.', dependencies=FALSE, build_vignettes=FALSE, upgrade = 'never')"
