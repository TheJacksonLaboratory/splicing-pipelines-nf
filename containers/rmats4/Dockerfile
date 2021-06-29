# continuumio/miniconda3:4.8.2
FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2

LABEL authors="phil@lifebit.ai laura.urbanski@jax.org" \
	description="Docker image containing rMATS v4.1.1"

COPY environment.yml /
RUN apt-get update && apt-get install -y procps && \
	conda env create -f /environment.yml && conda clean -a

# Make RUN commands use the new environment:
RUN echo "conda activate rmats4" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]	

ENV PATH /opt/conda/envs/rmats4/bin:$PATH
