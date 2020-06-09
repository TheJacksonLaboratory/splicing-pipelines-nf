FROM nfcore/base:1.9
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing dependencies to download reads from TCGA, GTEx and SRA"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/download_reads/bin:$PATH