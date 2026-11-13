FROM python:3.14

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends git wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Version floors match pyproject.toml. stsynphot 1.7.0 is not on PyPI
# yet; this URL is pull request 244. After the 1.7.0 tag, use
# "stsynphot>=1.7". See docs/getting_started.rst.
# NumPy, SciPy, and Astropy floors match that branch.
RUN pip install --no-cache-dir \
    "numpy>=2" "scipy>=1.14" matplotlib pandas "astropy>=7" \
    "synphot>=1.5" \
    "stsynphot @ git+https://github.com/wei-lingfeng/stsynphot_refactor.git@cache-grid-spectra"

RUN git clone https://github.com/MovingUniverseLab/spisea.git /SPISEA

ENV PYTHONPATH=/SPISEA

# Add data to cdbs folder
RUN mkdir /cdbs
WORKDIR /cdbs
RUN wget https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v10_sed.tar
RUN wget https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar
RUN tar -xvf hlsp_reference-atlases_hst_multi_everything_multi_v10_sed.tar && rm hlsp_reference-atlases_hst_multi_everything_multi_v10_sed.tar
RUN tar -xvf hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar && rm hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar

# Move data folders to cdbs
RUN mv /cdbs/grp/redcat/trds/comp /cdbs
RUN mv /cdbs/grp/redcat/trds/mtab /cdbs
RUN mv /cdbs/grp/redcat/trds/grid /cdbs

RUN mkdir /cdbs/models
WORKDIR /cdbs/models
RUN wget http://astro.berkeley.edu/~jlu/spisea/spisea_models.tar.gz && wget http://astro.berkeley.edu/~jlu/spisea/spisea_cdbs.tar.gz
RUN tar -xvf spisea_cdbs.tar.gz && tar -xvf spisea_models.tar.gz && rm spisea_cdbs.tar.gz && rm spisea_models.tar.gz

# stsynphot reads throughput/model grids from this path (same convention as legacy pysynphot).
ENV PYSYN_CDBS=/cdbs/models/cdbs/
ENV SPISEA_MODELS=/cdbs/models
