FROM continuumio/miniconda:4.6.14

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y

RUN conda install -c bioconda rmats=3.2.5 -y
RUN conda install -c anaconda pigz -y

ADD ./sampleCountsSave.sh /root/
RUN chmod 777 /root/sampleCountsSave.sh

ENV PATH="${PATH}:/root"
