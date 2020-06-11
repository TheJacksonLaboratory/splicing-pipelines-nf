# Running the pipeline

## For running the pipeline on Sumner
See [here](run_on_sumner.md)

## For running the pipeline on CloudOS
See [here](run_on_cloudos.md)

# Nextflow parameters

Nextflow parameters can be provided in one of two ways:
1) They can be specified in configuration files
2) They can be specified on the command-line

<details>
<summary>For example, all of the following are equivalent:</summary>

1) Config file
```groovy
params.reads      = '/path/to/reads.csv'
params.readlength = 48
params.singleEnd  = true
```

OR 

```groovy
params {
  reads      = '/path/to/reads.csv'
  readlength = 48
  singleEnd  = true
}
```

See [configuration scopes](https://www.nextflow.io/docs/latest/config.html#config-scopes) for more information on this^

2) Specifying parameters on the command-line
```bash
nextflow run main.nf --reads /path/to/reads.csv --readlength 48 --singleEnd true
```

</details>

Parameters specified on the command-line take precedence over those specified in configuration files. It is generally best-practice to have your parameters saved in a configuration file as this makes your analysis more reproducible if you need to run it again.

[Profiles](https://www.nextflow.io/docs/latest/en/latest/config.html#config-profiles) are configuration that can be included by specifying the profile name on the command-line. For example, `-profile sumner` to include configuration specific to JAX's HPC Sumner

## Types of Nextflow parameters

### String 

Strings can be specified using 'single' or "double quotes"

#### Paths

File paths can be any one of the following:
1) Local path - from directory that Nextflow is run in
2) Full path
3) Links - `https`, `ftp`, `s3` & `gs` links can all be used to specify input files provided that you have access to the file. Nextflow will automatically download these files into the `work` directory (in `work/stage`). On subsequent executions the pre-downloaded files will be used.

### Int

Integers can be specified without using quotes both in the configuration files and on the command-line

### Bool

Boolean parameters can be set to either `true` or `false`. Many of the parameters are initialised to `false` in [`nextflow.config`](../nextflow.config). You can set parameters to true on the command line just by using the flag. For example, just using `--singleEnd` will set the `singleEnd` parameter to true.

<details>
<summary>Side note:</summary>

However, be careful doing this as `--singleEnd false` will actually set the `singleEnd` parameter to the string `'false'` not the boolean `false`. Counterintuively, as this is a string that is present it actually mean that `singleEnd` will evaluate to true :satisfied:

This is another reason why it can be best to specify parameters in a confugration file rather than on the command-line
</details>

### Other

Memory units eg for `max_memory` can be specified in gigabytes eg `8.GB`

Time units eg for `max_time` can be specified in hours eg `2.h`

Both of these should be specified without quotes

## All available parameters
```
Main arguments:
  --reads                       Path to input data CSV file specifying the reads sample_id and path to FASTQ files (path)
  --gtf                         Path to GTF file (path)
  --star_index                  Path to STAR index (path)
  -profile                      Configuration profile to use. Can use multiple (comma separated, string)
                                Available: base, docker, sumner, test and more.

Reads:
  --rmats_pairs                 Path to file containing b1 & b2 samples names space seperated, one row for each rMATS comparison (path)
  --singleEnd                   Specifies that the input is single-end reads (bool)
  --stranded                    Specifies that the input is stranded (bool)
  --adapter                     Path to adapter file (path)
  --readlength                  Read length (int)
  --overhang                    Overhang (default = readlength - 1, int)
  --mismatch                    Mismatch (default = 2, int)

rMATS:
  --statoff                     Skip the statistical analysis (bool)
  --paired_stats                Use the paired stats model (bool)
  --novelSS                     Enable detection of novel splice sites (unannotated splice sites, bool)

Other:
  --assembly_name               Genome assembly name (available = 'GRCh38' or 'GRCm38', string)
  --test                        For running QC, trimming and STAR only (bool)
  --download_from               Database to download FASTQ/BAMs from (available = 'TCGA', 'GTEX' or 'SRA', string)
  --key_file                    For downloading reads, use TCGA authentication token (TCGA) or dbGAP repository key (GTEx, path)
  --max_cpus                    Maximum number of CPUs (int)
  --max_memory                  Maximum memory (memory unit)
  --max_time                    Maximum time (time unit)
  --skiprMATS                   Skip rMATS (bool)
  --skipMultiQC                 Skip MultiQC (bool)
  --outdir                      The output directory where the results will be saved (string)
```