FROM continuumio/miniconda3:23.3.1-0

LABEL maintainer="MarkerMAG"
LABEL description="MarkerMAG — Link MAGs to 16S rRNA genes (Nextflow workflow image)"

# Build with the repository root as the context:
#   docker build -f nextflow/Dockerfile -t markermag:latest .
WORKDIR /opt/markermag
COPY pyproject.toml setup.py MANIFEST.in LICENSE.txt /opt/markermag/
COPY bin/MarkerMAG /opt/markermag/bin/MarkerMAG
COPY MarkerMAG /opt/markermag/MarkerMAG
COPY nextflow/environment.yml /opt/markermag/nextflow/environment.yml

# environment.yml installs MarkerMAG-nf and the MATAM build/runtime dependencies.
WORKDIR /opt/markermag/nextflow
RUN conda env create --file environment.yml \
    && conda clean -afy

ENV PATH="/opt/conda/envs/markermag-nf/bin:${PATH}"

# Build the paired-read/Python 3.14 MATAM fork in the same Conda environment.
# Set MATAM_REF to a release tag or commit for reproducible production images.
ARG MATAM_REPOSITORY=https://github.com/xvazquezc/matam.git
ARG MATAM_REF=master
RUN git clone --recursive "${MATAM_REPOSITORY}" /opt/matam \
    && git -C /opt/matam checkout "${MATAM_REF}" \
    && git -C /opt/matam submodule update --init --recursive \
    && /opt/conda/envs/markermag-nf/bin/python /opt/matam/build.py

ENV PATH="/opt/matam/bin:/opt/conda/envs/markermag-nf/bin:${PATH}"

# ── NOTE: usearch is not freely redistributable.  ────────────────────────────
# Mount your licensed binary at runtime:
#   docker run -v /path/to/usearch:/usr/local/bin/usearch ...
# or pre-copy it into a derived image.

WORKDIR /data
