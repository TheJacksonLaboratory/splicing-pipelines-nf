# Running the pipeline

## For running the pipeline on Sumner
See [here](run_on_sumner.md)

## For running the pipeline on CloudOS
See [here](run_on_cloudos.md)

# Nextflow parameters

Nextflow parameters can be provided in one of two ways:
1) They can be specified in configuration files
2) They can be specified on the command-line

Parameters specified on the command-line take precedence over those specified in configuration files. It is generally best-practice to have your parameters saved in a configuration file as this makes your analysis more reproducible if you need to run it again.

[Profiles](https://www.nextflow.io/docs/latest/en/latest/config.html#config-profiles) are configuration that can be included by specifying the profile name on the command-line. For example, `-profile sumner` to include configuration specific to JAX's HPC Sumner

## Types of Nextflow parameters

### String 

String can be specified using single or double quotes

#### Paths

File paths can be any one of the following:
1) Local path - from directory that Nextflow is run in
2) Full path
3) Links - `https`, `ftp`, `s3` & `gs` links can all be used to specify input files provided that you have access to the file. Nextflow will automatically download these files into the `work` directory (in `work/stage`). On subsequent executions the pre-downloaded files will be used.

### Int

Integers can be specified without using quotes in the configuration files

### Bool

Boolean parameters can be set to either `true` or `false`. Many of the parameters are initialised to `false`. You can set parameters to true on the command line just by using the flag. For example, just using `--singleEnd` will set the `singleEnd` parameter to true.

### Other

Memory units eg for `max_memory` can be specified in gigabytes eg `8.GB`
Time units eg for `max_time` can be specified in hours eg `2.h`

## All available parameters
```
Main arguments:
  --reads                       Path to input data CSV file specifying the reads sample_id and path to FASTQ files
  --gtf                         Path to GTF file
  --star_index                  Path to STAR index
  -profile                      Configuration profile to use. Can use multiple (comma separated)
                                Available: docker, test and more.

Reads:
  --rmats_pairs                 Path to file containing b1 & b2 samples names space seperated, one row for each rMATS comparison
  --singleEnd                   Specifies that the input is single-end reads
  --stranded                    Specifies that the input is stranded
  --adapter                     Path to adapter file
  --readlength                  Read length (default = 48)
  --overhang                    Overhang (default = readlength - 1)
  --mismatch                    Mismatch (default = 2)

Other:
  --assembly_name               Genome assembly name
  --max_cpus                    Maximum number of CPUs
  --max_memory                  Maximum memory
  --max_time                    Maximum time
  --skiprMATS                   Skip rMATS
  --skipMultiQC                 Skip MultiQC
  --outdir                      The output directory where the results will be saved
```