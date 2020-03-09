FROM continuumio/miniconda:4.6.14

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y

RUN conda install -c bioconda sra-tools=2.9.6 -y
RUN conda install -c anaconda pigz -y
