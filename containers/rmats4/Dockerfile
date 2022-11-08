# continuumio/miniconda3:4.12.0
FROM continuumio/miniconda3@sha256:58b1c7df8d69655ffec017ede784a075e3c2e9feff0fc50ef65300fc75aa45ae

LABEL authors="phil@lifebit.ai laura.urbanski@jax.org" \
	description="Docker image containing rMATS v4.1.1"

RUN conda install -c conda-forge mamba -y

COPY environment.yml /
RUN apt-get update && apt-get install -y procps && \
	mamba env create -f /environment.yml && conda clean -a

# Make RUN commands use the new environment:
RUN echo "conda activate rmats4" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

ENV PATH /opt/conda/envs/rmats4/bin:$PATH
