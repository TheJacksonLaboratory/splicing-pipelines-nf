# Establishing setup before pipeline run

**Make sure you have the gcp-splice branch, it is the only branch tested for running the pipeline directly on GCP:**

### 1. JAVA version
Got to website  
```
https://www.oracle.com/java/technologies/downloads/#jdk19-mac
```
Download the correct DMG file based on the processor on your machine. After the dmg hasbeen downloaded install it into your machine. To check if the latest JDK version have been installed onto your machine. Type the command on your terminal:
```
java --version
```
The output should look somehting like:
```
java 17.0.5 2022-10-18 LTS
Java(TM) SE Runtime Environment (build 17.0.5+9-LTS-191)
Java HotSpot(TM) 64-Bit Server VM (build 17.0.5+9-LTS-191, mixed mode, sharing)
```

To make sure we have same javac and java version. Run the command:
```
javac --version
```
The output should look somehting like:
```
javac 17.0.5
```

### 2. Install google SDK on your local and make sure you have access to jax-cloudos-olgaanczukow project
It is good to follow the instruction provided at (https://cloud.google.com/sdk/docs/install-sdk) .... Once the installation is done
```
gcloud --version
gcloud init
gcloud auth login
gcloud auth application-default login
gcloud config configurations list
gcloud config set project jax-cloudos-olgaanczukow
```

### 3. Install miniconda
Please follow the instructions here (https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/macos.html)
In your terminal:
**Miniconda installer for macOS**
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
which conda
conda --version
```
### 4. Create and activate conda environment
```
conda env list
conda env create --name splicing-pipelines-nf -f containers/splicing-pipelines-nf/environment.yml
conda activate splicing-pipelines-nf
```
### 4. Install nextflow in the conda env
```
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
```
