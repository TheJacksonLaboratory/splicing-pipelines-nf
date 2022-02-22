FROM nfcore/base:1.9
LABEL authors="ines@lifebit.ai" \
      description="Docker image containing csvkit toolkit, including in2csv"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/csvkit/bin:$PATH