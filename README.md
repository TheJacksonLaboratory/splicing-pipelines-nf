# splicing-pipelines-nf
Repository for the Nextflow splicing pipeline

## Test the pipeline
Download the pipeline and test it on a minimal dataset with a single command:
```bash
nextflow run jacksonlabs/splicing-pipelines-nf -profile test,<docker/sumner>
```

## Example usage
Start running your own analysis:
```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome GRCh38 -profile <docker/sumner>
```