FROM continuumio/miniconda3@sha256:145567896379c33c89737aa1ab1389a71aead35b39f512118a13edcf9fe2a42c
LABEL authors="ines@lifebit.ai" \
      description="Docker image containing all requirements for the splicing-pipelines-nf"

ARG ENV_NAME="base"

RUN apt-get install -y procps

COPY environment.yml /
RUN conda env update -n ${ENV_NAME} -f environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ${ENV_NAME} > ${ENV_NAME}_exported.yml

# Initialise bash for conda
RUN conda init bash

COPY ./tagXSstrandedData.awk /usr/local/bin/
RUN chmod +x /usr/local/bin/tagXSstrandedData.awk

# Install the latest stringtie, gffread, gffcompare & deeptools
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.7.Linux_x86_64.tar.gz -O stringtie.tar.gz && \
    tar xvzf stringtie.tar.gz && mv stringtie-2.1.7.Linux_x86_64 stringtie && \
    rm stringtie.tar.gz && mv stringtie/prepDE.py stringtie/stringtie /usr/local/bin && \
    wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz -O gffread.tar.gz && \
    tar xvzf gffread.tar.gz && mv gffread-0.12.7.Linux_x86_64 gffread && \
    rm gffread.tar.gz && mv gffread/gffread /usr/local/bin/ && \
    wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz -O gffcompare.tar.gz && \
    tar xvzf gffcompare.tar.gz && mv gffcompare-0.12.6.Linux_x86_64 gffcompare && \
    rm gffcompare.tar.gz && mv gffcompare/gffcompare gffcompare/trmap /usr/local/bin/

# Install gawk to fix https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/issues/120
RUN apt-get update && apt-get install -y gawk