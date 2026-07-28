FROM continuumio/miniconda3:23.3.1-0

LABEL maintainer="MarkerMAG"
LABEL description="MarkerMAG — Link MAGs to 16S rRNA genes (Nextflow workflow image)"

# Build with the repository root as the context:
#   docker build -f nextflow/Dockerfile -t markermag:latest .
WORKDIR /opt/markermag
COPY setup.py MANIFEST.in LICENSE.txt /opt/markermag/
COPY bin/MarkerMAG /opt/markermag/bin/MarkerMAG
COPY MarkerMAG /opt/markermag/MarkerMAG
COPY nextflow/environment.yml /opt/markermag/nextflow/environment.yml

# environment.yml installs the local package in editable mode from "..".
WORKDIR /opt/markermag/nextflow
RUN conda env update -n base --file environment.yml \
    && conda clean -afy

# ── NOTE: usearch is not freely redistributable.  ────────────────────────────
# Mount your licensed binary at runtime:
#   docker run -v /path/to/usearch:/usr/local/bin/usearch ...
# or pre-copy it into a derived image.

WORKDIR /data
