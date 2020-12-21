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
Input files:
  --reads                       Path to reads.csv file, which specifies the sample_id and path to FASTQ files for each read or read pair (path).
                                This file is used if starting at beginning of pipeline. 
                                (default: no reads.csv)
  --bams                        Path to bams.csv file which specifies sample_id and path to BAM and BAM.bai files (path)
                                This file is used if starting pipeline at Stringtie.
                                (default: no bams.csv)
  --rmats_pairs                 Path to rmats_pairs.txt file containing b1 (and b2) samples names (path)
                                (default: no rmats_pairs specified) 
  --run_name                    User specified name used as prefix for output files
                                (defaut: no prefix, only date and time)
  --download_from               Database to download FASTQ/BAMs from (available = 'TCGA', 'GTEX' or 'GEN3-DRS', 'SRA') (string)
                                (default: false)
  --key_file                    For downloading reads, use TCGA authentication token (TCGA) or dbGAP repository key (GTEx, path)
  				or credentials.josn file in case of 'GEN3-DRS'
                                (default: false)
                                
Main arguments:
  --gtf                         Path to reference GTF file (path)
                                (default: no gtf specified) 
  --assembly_name               Genome assembly name (available = 'GRCh38' or 'GRCm38', string)
                                (default: false)
  --star_index                  Path to STAR index (path)
                                (default: read length)
  --singleEnd                   Specifies that the input is single-end reads (bool)
                                (default: false)
  --stranded                    Specifies that the input is stranded ('first-strand', 'second-strand', false (aka unstranded))
                                (default: 'first-strand')
  -profile                      Configuration profile to use. Can use multiple (comma separated, string)
                                Available: base, docker, sumner, test and more.
  --readlength                  Read length - Note that all reads will be cropped to this length(int)
                                (default: no read length specified)
                                
Trimmomatic: 
  --minlen                      Drop the read if it is below a specified length (int)
				Default parameters turn on --variable-readlength
				To crop all reads and turn off, set minlen = readlength (NOTE: this will turn off soft clipping)                                
                                (default: 20)
  --slidingwindow               Perform a sliding window trimming approach (bool)
                                (default: true)
  --adapter                     Path to adapter file (path)  
                                (default: TruSeq3 for either PE or SE, see singleEnd parameter)
                                
Star:                    
  --mismatch                    Number of allowed mismatches per read (SE) or combined read (PE) (int)
                                SE ex. read length of 50, allow 2 mismatches per 50 bp
                                PE ex. read length of 50, allow 2 mismatches per 100 bp 
                                (default: 2)
  --overhang                    Overhang (int)
                                (default: readlength - 1)
  --filterScore                 Controls --outFilterScoreMinOverLread and outFilterMatchNminOverLread
                                (default: 0.66)
  --sjdbOverhangMin             Controls --alignSJDBoverhangMin (int)
                                (default: 3)
  --star_memory                 Max memory to be used by STAR to sort BAM files.
                                (default: Available task memory)

rMATS:                              
  --statoff                     Skip the statistical analysis (bool)
                                If using only b1 as input, this must be turned on.
                                (default: false)
  --paired_stats                Use the paired stats model (bool)
                                (default: false)
  --novelSS                     Enable detection of unnanotated splice sites (bool)
                                (default: false)
  --mil                         Minimum Intron Length. Only impacts --novelSS behavior (int)
                                (default: 50)
  --mel                         Maximum Exon Length. Only impacts --novelSS behavior (int)
                                (default: 500)

Other:
  --test                        For running trim test (bool)
                                (default: false)
  --max_cpus                    Maximum number of CPUs (int)
                                (default: 72)  
  --max_memory                  Maximum memory (memory unit)
                                (default: 760)
  --max_time                    Maximum time (time unit)
                                (default: 72.h)
  --skiprMATS                   Skip rMATS (bool)
                                (default: false)
  --skipMultiQC                 Skip MultiQC (bool)
                                (default: false)
  --outdir                      The output directory where the results will be saved (string)
                                (default: directory where you submit the job)
  --gc_disk_size                Only specific to google-cloud executor. Adds disk-space for few aggregative processes.
                                (default: "200 GB" based on 100 samples. Simply add 2 x Number of Samples)
  --mega_time                   Sets time limit for processes withLabel 'mega_memory' in the main.nf using the base.config (time unit)
                                Make sure '#SBATCH -t' in 'main.pbs' is appropriately set if you are changing this parameter.
                                (default: 20.h)

	    
```

## Run with data from AnviL Gen3-DRS

You will be needing two things from - https://gen3.theanvil.io/

1. manifest file
2. credentials file

Original downloaded `manifest.json` file need to be converted into `manifest.csv` in order to be accepted in `--reads`, for doing that you can do this - 

```bash
pip install csvkit
in2csv manifest.json > manifest.csv
```

NOTE: Make sure the `manifest.csv` file have five columns, Check from [examples](../examples/gen3/)

Downloaded `credentials.json` file can be provided in `--key` param.

NOTE: Make sure `credentials.json` is a latest one. They have expiry dates when you download.

If you running with AnviL Gen3-DRS files you also need to provide a Genome fasta file with `--genome_fasta`, which will be used to convert CRAM files to BAM format.

For a minimal params list check [gen3_drs.config](../conf/examples/GEN3_DRS_config.md)

### Extract based on a bam query list

If you have a list of bam file names of interest, extract the manifest file - 

```bash
# Get all the bam files name into a txt file
cut -d, -f4 query_List.csv > query_list.txt
# Extract those bam files list from manifest.csv
grep -f query_list.txt -i manifest.csv > manifest.csv
```

