FROM registry.codeocean.com/codeocean/r-base:3.5.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      "cmake=3.10.2-1ubuntu2" \
      "libarmadillo-dev=1:8.400.0+dfsg-2" \
      "libboost-all-dev=1.65.1.0ubuntu1" \
      "liblapack-dev=3.7.1-4ubuntu1" \
      "libopenblas-dev=0.2.20+ds-4" \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(warn=2); devtools::install_version("ROCR", version = "1.0-7")'
RUN Rscript -e 'options(warn=2); devtools::install_version("circlize", version = "0.4.6")'
RUN Rscript -e 'options(warn=2); devtools::install_version("dendsort", version = "0.3.3")'
RUN Rscript -e 'options(warn=2); devtools::install_version("glmnet", version = "2.0-18")'
RUN Rscript -e 'options(warn=2); devtools::install_version("gridExtra", version = "2.3")'
RUN Rscript -e 'options(warn=2); devtools::install_version("openxlsx", version = "4.1.0.1")'
RUN Rscript -e 'options(warn=2); devtools::install_version("reshape", version = "0.8.8")'
RUN Rscript -e 'options(warn=2); devtools::install_version("wordcloud", version = "2.6")'

RUN Rscript -e ' \
      source("http://bioconductor.org/biocLite.R"); \
      biocLite(c( \
        "ComplexHeatmap" \
      ), suppressUpdates = TRUE)'
