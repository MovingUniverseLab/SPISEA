FROM python:3.14

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends git wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Version floors match pyproject.toml.
RUN pip install --no-cache-dir \
    "numpy>=2" "scipy>=1.14" matplotlib pandas "astropy>=7" \
    "synphot>=1.9" "stsynphot>=1.7"

RUN git clone https://github.com/MovingUniverseLab/spisea.git /SPISEA

ENV PYTHONPATH=/SPISEA

# Evolution grids and SPISEA R=250 atmosphere catalogs (grid_version 3.1).
# spisea_models.tar.gz unpacks to evolution/.
# spisea_cdbs.tar.gz unpacks to cdbs/grid/.
# Peak disk during this step is about 30 GB (tarball plus unpacked tree).
RUN mkdir -p /models
WORKDIR /models
RUN wget --progress=dot:giga \
        https://w.astro.berkeley.edu/~jlu/spisea/spisea_models.tar.gz \
        https://w.astro.berkeley.edu/~jlu/spisea/spisea_cdbs.tar.gz \
    && tar --warning=no-unknown-keyword -xzf spisea_models.tar.gz \
    && tar --warning=no-unknown-keyword -xzf spisea_cdbs.tar.gz \
    && rm spisea_models.tar.gz spisea_cdbs.tar.gz

# STScI Reference Atlases. Filter throughputs (comp/, mtab/) and the Phoenix
# grid have to sit in the same cdbs tree as the SPISEA catalogs.
# synphot1_throughput-master.tar is the v20 atlas (Roman filters, 2024-03-22+).
# synphot5_phoenix-models.tar is the Phoenix grid.
RUN mkdir -p /tmp/synphot1 /tmp/synphot5 \
    && wget --progress=dot:giga -O /tmp/synphot1.tar \
        https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v20_sed.tar \
    && wget --progress=dot:giga -O /tmp/synphot5.tar \
        https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_phoenix-models_multi_v3_synphot5.tar \
    && tar -xf /tmp/synphot1.tar -C /tmp/synphot1 \
    && tar -xf /tmp/synphot5.tar -C /tmp/synphot5 \
    && mv /tmp/synphot1/grp/redcat/trds/comp /models/cdbs/comp \
    && mv /tmp/synphot1/grp/redcat/trds/mtab /models/cdbs/mtab \
    && cp -a /tmp/synphot5/grp/redcat/trds/grid/. /models/cdbs/grid/ \
    && rm -rf /tmp/synphot1 /tmp/synphot5 /tmp/synphot1.tar /tmp/synphot5.tar

# stsynphot reads throughput and atmosphere grids from PYSYN_CDBS.
# SPISEA_MODELS is the directory that contains evolution/.
ENV PYSYN_CDBS=/models/cdbs
ENV SPISEA_MODELS=/models
