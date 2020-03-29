# splicing-pipelines-nf
Anczukow-lab repository for the Nextflow splicing pipeline with `rMATS 4.0.2`.

## Sumner (JAX HPC) execution

Check out the source code from `github`.

```bash
git clone https://github.com/TheJacksonLaboratory/splicing-pipelines-nf.git
cd splicing-pipelines-nf
```

Then to execute on sumner, a user does the following bash script with the `slurm` `sbatch` command:

```bash
sbatch main.pbs
```
Progress can be checked using:

```bash
squeue | grep [username]
```
or by tailing the output file:

```bash
tail -f splicing.[jobnumber].out
```

## Diving into the details of the Example execution

The `main.pbs` executes an analysis comparing `MYC_MCF10A_0h` with 3 replicates and `MYC_MCF10_8h`.
The details of what needs to be configured to do this comparison analysis are found in three files:

All the analyses are kept in `analyses` subdirectory (by convention).   Encoded in the file name is the metadata details that outline the comparison that is being completed.  In this case capturing the statement above (`MYC_MCF10A_0h` vs `MYC_MCF10A_8h`).

These files can be specified via command line or via a config file.

To specify via command line:

* `--reads analysis/MYC_MCF10A_0h_vs_MYC_MCF10A_8h.csv`
    This file contains the `sample_id` a short name uniquely defines the sample within this comparsion
    comma seperated with the complete path for the `left` and `right` `fastqs`.   
    
* `--b1 analysis/MYC_MCF10A_0h_vs_MYC_MCF10A_8h_b1.txt`
    This is a comma separated file containing 1 to many replicates for the `case` in the example.
    
* `--b2 analysis/MYC_MCF10A_0h_vs_MYC_MCF10A_8h_b2.txt`
    This is a comma separated file containing 1 to many replicates for the `control` in the example.

## NextFlow Config Order

For the splicing pipeline we have:

### nextflow.config

* `nextflow.config` is the main config file found at the `root directory`. It includes defaults for all parameters, containers & links to other configs/profiles
* `base.config` to specify resources
* `igenomes.config` specifying the path to the igenomes resources available
    To mirror the setup on the CloudOS Universal Research Platform, adopting this standard setup.
    * Currently in `/projects/adeslat/igenomes/` `GRCh38` has been downloaded
    * `STARIndex` downloaded from the `AWS S3 Budget` Amazon has donated to host these standards
* `igenomes.config` to specify path to files for a given reference genome

## All parameters
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
