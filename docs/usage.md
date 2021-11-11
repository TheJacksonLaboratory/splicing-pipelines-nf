# Running the pipeline

## For running the pipeline on Sumner
See [here](run_on_sumner.md)

## For running the pipeline on CloudOS
See [here](run_on_cloudos.md)

## For running the pipeline locally
See [here](run_locally.md)

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
    --reads                     Path to reads.csv file, which specifies the sample_id and path to FASTQ files
                                for each read or read pair (path).
                                This file is used if starting at beginning of pipeline. It can be file paths,
                                s3 links or ftp link.
                                (default: no reads.csv)
    --bams                      Path to bams.csv file which specifies sample_id and path to BAM and BAM.bai
                                files (path)
                                If this file is provided, pipeline will start at Stringtie (and proceed through
                                rMATS and post processing).
                                (default: no bams.csv)
    --rmats_pairs               Path to rmats_pairs.txt file containing b1 (and b2) samples names (path)
                                (default: no rmats_pairs specified)
    --run_name                  User specified name used as prefix for output files
                                (defaut: no prefix, only date and time)
    --download_from             Database to download FASTQ/BAMs from (available = 'TCGA', 'GTEX' or 'GEN3-DRS',
                                'SRA', 'FTP') (string)
                                false should be used to run local files on the HPC (Sumner).
                                'TCGA' can also be used to download GDC data including HCMI data.
                                (default: false)
    --key_file                  For downloading reads, use TCGA authentication token (TCGA) or dbGAP repository
                                key (GTEx, path) or credentials.json file in case of 'GEN3-DRS'
                                (default: false)

Main arguments:
    --gtf                       Path to reference GTF file (path)
                                (default: no gtf specified)
    --assembly_name             Genome assembly name (available = 'GRCh38' or 'GRCm38', string)
                                (default: false)
    --star_index                Path to STAR index (path)
                                Star indices must be generated prior to run (with correct STAR version)
                                (default: false)
    --singleEnd                 Specifies that the input is single-end reads (bool)
                                This parameter also automatically establishes the path to the SE or PE adapters.
                                For PE, set to false.
                                (default: false)
    --stranded                  Specifies that the input is stranded ('first-strand', 'second-strand',
                                false (aka unstranded))
                                'first-strand' refers to RF/fr-firststrand in this pipeline.
                                (default: 'first-strand')
    --readlength                Read length - Note that all reads will be cropped to this length(int)
                                (default: no read length specified)
    -profile                      Configuration profile to use. Can use multiple (comma separated, string)
                                On sumner, this should be set in the main.pbs or as a command-line parameter.
                                Profile can only be activated from the command line.
                                Available: base, docker, sumner, test and more.

Trimmomatic:
    --minlen                    Drop the read if it is below a specified length (int)
                                Default parameters turn on --variable-readlength
                                To crop all reads and turn off --variable-readlength, set minlen = readlength
                                (default: 20)
    --slidingwindow             Perform a sliding window trimming approach (bool)
                                (default: true)
    --adapter                   Path to adapter file (path)
                                (default: TruSeq3 for either PE or SE, see singleEnd parameter)

Star:
    --mismatch                  Number of allowed mismatches per read (SE) or combined read (PE) (int)
                                SE ex. read length of 50, allow 2 mismatches per 50 bp
                                PE ex. read length of 50, allow 2 mismatches per 100 bp
                                (default: 2)
    --overhang                  Overhang (int)
                                (default: readlength - 1)
    --filterScore               Controls --outFilterScoreMinOverLread and outFilterMatchNminOverLread
                                For TCGA values:
                                https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
                                (default: 0.66)
    --sjdbOverhangMin           Controls --alignSJDBoverhangMin (int)
                                For TCGA values:
                                https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
                                (default: 3)
    --soft_clipping             Enables soft clipping (bool)
                                If true, the STAR parameter will be --alignEndsType 'Local' and the rMATS parameter
                                --allow-clipping will be added.
                                If false, the STAR parameter will be --alignEndsType 'EndToEnd' and no rMATS
                                parameter is added.
                                NOTE: Soft Clipping will cause read lengths to be variable, so turn soft_clipping
                                off if reads need to be same length. Variable read length parameter is turned on
                                in rMATS when minlen does not equal readlength.
                                (default: true)
    --save_unmapped             Save unmapped and partially mapped reads in separate file (bool)
                                (default: false)
    --star_memory               Max memory to be used by STAR to sort BAM files.
                                (default: Available task memory)

rMATS:
    --statoff                   Skip the statistical analysis (bool)
                                If using only b1 as input, this must be turned on.
                                (default: false)
    --paired_stats              Use the paired stats model (bool)
                                (default: false)
    --novelSS                   Enable detection of unnanotated splice sites (bool)
                                (default: false)
    --mil                       Minimum Intron Length. Only impacts --novelSS behavior (int)
                                (default: 50)
    --mel                       Maximum Exon Length. Only impacts --novelSS behavior (int)
                                (default: 500)

Other:
    --test                      For running trim test (bool)
                                To run the first half of the pipeline (through STAR), set test = true.
                                (default: false)
    --max_cpus                  Maximum number of CPUs (int)
                                (default: 72)
    --max_memory                Maximum memory (memory unit)
                                (default: 760.GB)
    --max_time                  Maximum time (time unit)
                                (default: 72.h)
    --skiprMATS                 Skip rMATS (bool)
                                (default: false)
    --skipMultiQC               Skip MultiQC (bool)
                                (default: false)
    --outdir                    The output directory where the results will be saved (string)
                                On Sumner, this must be set in the main.pbs or via command line.
                                NF_splicing_pipeline.config will not overwrite main.pbs.
                                (default: <directory where you submit the job>/results)
    --mega_time                 Sets time limit for processes withLabel 'mega_memory' in the main.nf using the
                                base.config (time unit)
                                Make sure '#SBATCH -t' in 'main.pbs' is appropriately set if you are changing this parameter.
                                (default: 20.h)
    --gc_disk_size              Only specific to google-cloud executor. Adds disk-space for few aggregative processes.
                                (default: "200 GB" based on 100 samples. Simply add 2 x Number of Samples)
    --debug                     This option will enable echo of script execution into STDOUT with some additional
                                resource information (such as machine type, memory, cpu and disk space)
                                (default: false)
    --error_strategy            Mode of pipeline handling failed processes.
                                Possible values: 'terminate', 'finish', 'ignore', 'retry'.
                                Check nextflow documnetation for detailed descriptions of each mode:
                                https://www.nextflow.io/docs/latest/process.html#process-page-error-strategy
                                Set this parameter in the main.pbs, on the command line, or see NF_splicing_pipeline.config
                                example (does not work like normal config param)
                                This does not overwrited CloudOS config, which is set to:
                                'errorStrategy = { task.exitStatus in [3,9,10,14,143,137,104,134,139] ? 'retry': 'ignore'}
                                (default (non-cloudos): 'finish')
    --cleanup                   This option will enable nextflow work folder cleanup upon pipeline successfull
                                completion. All intermediate files from nexftlow processes' workdirs will be
                                cleared, staging folder with staged files will not be cleared.
                                If pipeline is completed with errors or interrupted cleanup will not be executed.
                                Following successfull run resumed from the failed run with --cleanup option enabled
                                will only clear folders of processess created in the latest run, it will not clear
                                cached folders coming from previous pipleine runs.
                                Set this parameter in the main.pbs, on the command line, or see NF_splicing_pipeline.config
                                example (does not work like normal config param)
                                (default non-cloudos: true; cloudos: false)
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
cut -d, -f4 query_list.csv > bam_files_list.txt
# Extract those bam files list from manifest.csv
grep -f bam_files_list.txt -i manifest.csv > manifest.csv
```

Here `query_list.csv` should look something like - 

```csv
file_name,sequencing_assay,data_format,file_name,sample_id,participant_id,tissue,age,gender
GTEX-PPPP-XXX-XX-XXXXX,RNA-Seq,bam,GTEX-PPPP-XXX-XX-XXXXX.Aligned.sortedByCoord.out.patched.md.bam,GTEX-PPPP-XXX-XX-XXXXX,GTEX-PPPP,Breast,21,Female
```
