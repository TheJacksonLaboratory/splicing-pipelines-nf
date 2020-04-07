#!/bin/bash
#
# Pre-requisite - these data need to be downloaded by an authorized and authenticated user 
# on the google cloud platform
#
# gcloud auth login
# then this command may be run.
#

gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039513_1.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039513_2.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039509_1.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039509_2.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039517_1.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039517_2.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039508_1.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039508_2.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039512_1.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039512_2.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039516_1.fastq.gz .
gsutil -m cp gs://cloudosinputdata/inputs/splicing-pipelines-nf/samples/SRR1039516_2.fastq.gz .

