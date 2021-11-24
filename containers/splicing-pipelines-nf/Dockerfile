FROM nfcore/base:1.9
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing all requirements for the splicing-pipelines-nf"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/splicing-pipelines-nf/bin:$PATH

COPY ./tagXSstrandedData.awk /usr/local/bin/
RUN chmod +x /usr/local/bin/tagXSstrandedData.awk

# Install the latest stringtie, gffread & gffcompare
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
RUN apt-get install -y gawk