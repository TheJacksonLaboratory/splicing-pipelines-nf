# Run locally

This guide provides information needed to run the pipeline in a local linux environment.

## 0) Install Nextflow, Docker and/or Singularity

### 0.1) Install Nextflow
```bash
wget -qO- https://get.nextflow.io | bash && sudo mv nextflow /usr/bin/
```

### 0.1) Install Docker
See https://docs.docker.com/engine/install/ubuntu/


### 0.2) Install Singularity
See https://sylabs.io/guides/3.8/user-guide/quick_start.html and https://sylabs.io/guides/3.0/user-guide/installation.html

Commands to install on CentOS:
```bash
sudo yum update -y && \
sudo yum groupinstall -y 'Development Tools' && \
sudo yum install -y \
    openssl-devel \
    libuuid-devel \
    libseccomp-devel \
    wget \
    squashfs-tools

export VERSION=1.16.6 OS=linux ARCH=amd64 && # adjust this as necessary \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

go get -u github.com/golang/dep/cmd/dep


export VERSION=3.8.0 && # adjust this as necessary \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
    tar -xzf singularity-ce-${VERSION}.tar.gz && \
    cd singularity-ce-${VERSION}

./mconfig && \
    make -C builddir && \
    sudo make -C builddir install

singularity help

```

## 1) Clone the git repository to local machine

```bash
git clone https://github.com/TheJacksonLaboratory/splicing-pipelines-nf
cd splicing-pipelines-nf

```

## 2) Run the pipeline

### 2.1) Quick test with Docker container engine
```bash
nextflow run . -profile ultra_quick_test,docker --cleanup
```

```bash
nextflow run . -profile docker --singleEnd=true --reads="examples/testdata/single_end/tiny_reads.csv" --gtf='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/genes.gtf' --star_index='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/star_2.7.9a_yeast_chr_I.tar.gz' --test=true --readlength=48 --overhang=100 --max_cpus=2 --max_memory='6.GB' --max_time='48.h'
```

### 2.2) Quick test with Singularity container engine

```bash
nextflow run . -profile ultra_quick_test,singularity --cleanup
```

```bash
nextflow run . -profile singularity --singleEnd=true --reads="examples/testdata/single_end/tiny_reads.csv" --gtf='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/genes.gtf' --star_index='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/star_2.7.9a_yeast_chr_I.tar.gz' --test=true --readlength=48 --overhang=100 --max_cpus=2 --max_memory='6.GB' --max_time='48.h'
```

### 2.3) SRA quick test with Docker container engine
```bash
nextflow run . -profile sra_test,docker --cleanup
```

```bash
nextflow run . -profile docker --singleEnd=true --reads="examples/testdata/sra/sra.csv" --download_from='FTP' --gtf='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/genes.gtf' --star_index='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/star_2.7.9a_yeast_chr_I.tar.gz' --test=true --readlength=500 --overhang=100 --max_cpus=2 --max_memory='6.GB' --max_time='48.h'
```

### 2.4) SRA quick test with Singularity container engine

```bash
nextflow run . -profile sra_test,singularity --cleanup
```

```bash
nextflow run . -profile docker --singleEnd=true --reads="examples/testdata/sra/sra.csv" --download_from='FTP' --gtf='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/genes.gtf' --star_index='https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/star_2.7.9a_yeast_chr_I.tar.gz' --test=true --readlength=500 --overhang=100 --max_cpus=2 --max_memory='6.GB' --max_time='48.h'
```