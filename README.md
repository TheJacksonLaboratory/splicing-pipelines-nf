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
nextflow run main.nf --reads examples/splicing_pipeline_testdata/single_end_reads.csv --genome GRCh38 -profile <docker/sumner>
```

## All parameters
```
Mandatory arguments:
  --reads                       Path to input data CSV file specifying the reads sample_id and path to FASTQ files
  --singleEnd                   Specifies that the input is single-end reads
  --genome                      Name of iGenomes reference
  -profile                      Configuration profile to use. Can use multiple (comma separated)
                                Available: docker, test and more.

References:
  --gtf                         Path to GTF file
  --star_index                  Path to STAR index

Reads:
  --stranded                    Specifies that the input is stranded
  --adapter                     Path to adapter file
  --readlength                  Read length (default = 48)
  --overhang                    Overhang (default = readlength - 1)
  --mismatch                    Mismatch (default = 2)

Other:
  --max_cpus                    Maximum number of CPUs
  --max_memory                  Maximum memory
  --max_time                    Maximum time
  --outdir                      The output directory where the results will be saved
  --igenomes_base               Path to iGenomes base directory (default = 's3://ngi-igenomes/igenomes/')
```