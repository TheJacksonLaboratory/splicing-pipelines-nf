FROM nfcore/base:1.9
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing all requirements for the splicing-pipelines-nf"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/splicing-pipelines-nf/bin:$PATH