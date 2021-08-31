# continuumio/miniconda3:4.8.2
FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2

LABEL authors="phil@lifebit.ai" \
      description="Docker image containing dependencies to download reads from TCGA, GTEx and SRA"

RUN apt-get install procps -y

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/download_reads/bin:$PATH