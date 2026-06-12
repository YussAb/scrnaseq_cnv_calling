FROM rocker/r-ver:4.3.3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        cmake \
        gfortran \
        git \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgsl-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff5-dev \
        libxml2-dev \
        make \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("remotes", "parallelDist", "dlm", "gplots", "RColorBrewer", "mixtools", "MCMCpack"), repos = "https://cloud.r-project.org")'

RUN Rscript -e 'remotes::install_github("navinlabcode/copykat", upgrade = "never")'

RUN Rscript -e 'library(copykat); data("full.anno", package = "copykat"); data("DNA.hg20", package = "copykat"); data("cyclegenes", package = "copykat"); test_gene <- full.anno$hgnc_symbol[which(!is.na(full.anno$hgnc_symbol))[1]]; test_mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list(test_gene, "cell_1")); stopifnot(nrow(copykat::annotateGenes.hg20(test_mat, ID.type = "S")) == 1); cat("copykat", as.character(utils::packageVersion("copykat")), "\n")'
