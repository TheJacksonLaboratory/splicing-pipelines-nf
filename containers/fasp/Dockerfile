# continuumio/miniconda3:4.8.2
FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2

################## METADATA ######################
LABEL software="ga4gh/fasp" \
      about.home="https://github.com/lifebit-ai/fasp-scripts" \
      maintainer="Sangram Keshari Sahu <sangram@lifebit.ai>"

################## INSTALLATION ######################

COPY requirements.txt /

RUN apt-get update \
  && apt-get install procps -y \
  && pip install -r requirements.txt \
  && git clone https://github.com/lifebit-ai/fasp-scripts --branch gen3-drs-gtex-download \ 
  && cd fasp-scripts \
  && python setup.py install \
  && chmod +x fasp/scripts/* \
  && conda install samtools=1.11 -c bioconda -c conda-forge \
  && conda clean -a
  
ENV PATH /fasp-scripts/fasp/scripts:$PATH
