## Based on the RStudio image with R 4.6 (matches rproject.toml) ----

FROM rocker/rstudio:4.6.0

LABEL org.opencontainers.image.authors="Chung Shing Rex Ha <hachungshingrex@gmail.com>"


## Install system dependencies (from `rv sysdeps`) ----

RUN apt-get update -yq \
 && apt-get install --no-install-recommends -yq \
    cmake curl libcurl4-openssl-dev libfontconfig1-dev libfreetype6-dev \
    libfribidi-dev libgit2-dev libglpk-dev libharfbuzz-dev libicu-dev \
    libpng-dev libtiff-dev libuv1-dev libwebp-dev libx11-dev libxml2-dev \
    zlib1g-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


## Install rv (pinned to the version used locally) ----

ARG RV_VERSION=0.23.1
RUN curl -fsSL https://github.com/A2-ai/rv/releases/download/v${RV_VERSION}/rv-v${RV_VERSION}-$(uname -m)-unknown-linux-gnu.tar.gz \
    | tar -xz -C /usr/local/bin


## Install R packages (cached until the lockfile changes) ----

ENV FOLDER=/home/rstudio/
WORKDIR $FOLDER

COPY --chown=rstudio:rstudio rproject.toml rv.lock $FOLDER
USER rstudio
RUN rv sync && rv activate
USER root


## Copy local project ----

COPY --chown=rstudio:rstudio . $FOLDER
