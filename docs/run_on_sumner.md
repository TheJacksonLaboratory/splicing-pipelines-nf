## Steps needed to run a NextFlow workflow on JAX internal HPC Sumner
1) Login to sumner
```bash
ssh adeslat@login.sumner.jax.org
```
2) Change to a directory where you can checkout the code on github
```bash
cd /projects/adeslat
```
3) Clone the repository required
```bash
git clone https://github.com/TheJacksonLaboratory/splicing-pipelines-nf.git
```
4) Change into the cloned repository
```bash
cd splicing-pipelines-nf
```
5) Download or use your already downloaded nextflow process
```bash
curl -fsSL get.nextflow.io | bash
```
6) Execute the nextflow pipeline.
```bash
./nextflow run main.nf -profile test,sumner 
```

## Troubleshooting
1) Log into an interactive node
```bash
srun -n 1 --mem 1000 --pty /bin/bash
```
2) Load the singularity module
```bash
module load singularity
```
