FROM continuumio/miniconda3@sha256:a2e6aa4cd0b6dd696ae9e3e5732943250a977ab3a42b2fe5fb7ef0c19d2d9f16
LABEL authors="ines@lifebit.ai" \
      description="Docker image containing all requirements for the splicing-pipelines-nf"

ARG ENV_NAME="base"

RUN apk add --no-cache bash=5.0.17-r0 procps=3.3.16-r0


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
RUN apk add --no-cache gawk