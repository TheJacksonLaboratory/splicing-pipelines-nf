# Running the pipeline

## For running the pipeline on Sumner
See [here](run_on_sumner.md)

## For running the pipeline on CloudOS
See [here](run_on_cloudos.md)

## All available parameters
```
Main arguments:
  --reads                       Path to input data CSV file specifying the reads sample_id and path to FASTQ files
  --singleEnd                   Specifies that the input is single-end reads
  --genome                      Name of iGenomes reference
  -profile                      Configuration profile to use. Can use multiple (comma separated)
                                Available: docker, test and more.

References:
  --gtf                         Path to GTF file
  --star_index                  Path to STAR index
  --gencode_gtf                 Path to gencode file (annotation gtf)

Reads:
  --stranded                    Specifies that the input is stranded
  --b1                          Path to rMATs b1 file containing sample names
  --b2                          Path to rMATs b2 file containing sample names
  --adapter                     Path to adapter file
  --readlength                  Read length (default = 48)
  --overhang                    Overhang (default = readlength - 1)
  --mismatch                    Mismatch (default = 2)

Other:
  --max_cpus                    Maximum number of CPUs
  --max_memory                  Maximum memory
  --max_time                    Maximum time
  --skiprMATS                   Skip rMATS
  --skipMultiQC                 Skip MultiQC
  --outdir                      The output directory where the results will be saved
  --igenomes_base               Path to iGenomes base directory (default = 's3://ngi-igenomes/igenomes/')
```