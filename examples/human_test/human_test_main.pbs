#!/bin/bash
#SBATCH -o splicing.%j.out
#SBATCH -e splicing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$USER@jax.org
#SBATCH --mem=20000
#SBATCH --cpus-per-task=4
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 2-09:00:00

### to run this test: sbatch /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/examples/human_test/human_test_main.pbs
 
cd $SLURM_SUBMIT_DIR
date;hostname;pwd

module load singularity

export NXF_VER=20.04.1

curl -fsSL get.nextflow.io | bash

./nextflow run /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/main.nf \
	-config /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/conf/examples/human_test.config  \
	--outdir ${SLURM_SUBMIT_DIR} \
	-profile sumner -resume \
	-with-report human_test.html \
	-with-timeline human_test_timeline.html 
