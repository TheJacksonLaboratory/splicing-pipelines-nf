# Changelog

## Current Pipeline improvements
#### Improvements:
 - Removes network and subnetwork parameters because CloudOS now sets this automatically (#315).

#### Fixes:

#### Updates:

#### Maintenance:

## v2.1 - Pipeline improvements

#### Improvements:
- removed fasta file requirement from `Gen3-DRS` option (#297)
- The `google-lifesciences` has been defined as the executor when google profile is used. (#305)
- `GTEX`downlod option has been removed due to being obsolete (#299)
- `Gen3-DRS` was renamed to `GTEX` as it's now the only way to download GTEX file (#299)
- Added a new `--manifest` parameter with the input being the .json manifest file downloaded from GTEX (#304)
- Added `--reads` parameter which receives a CSV file with samples for which the analysis will be limited to when using `GTEX` (#304)


#### Fixes:
- Fix pipeline crach when using `Gen3-DRS` input due to `env(...)`format being used in process output (#294)
- Changed `gtex`input to support new manifest file format (#296)

#### Updates:
- None performed

#### Maintenance:
- None added


## v2.0 - Pipeline improvements
#### Improvements:
 - Adds saving of all the process .command* log files to results/process-logs folder (#251)
 - Adds pipeline workdir `--cleanup` option to clear all intermediate files on pipeline successful completion (true by default, false for CloudOS) (#238, #284, [089d6e3](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/pull/245/commits/3b71e038b186bb2bc92debacb02aede7b5dae917))
 - Adds pipeline `--error_strategy` parameter to be able to specify pipeline error strategy directly from command line (doesn't work if specified in config linked by `-c` or `-config` nextflow params) (#267)
 - Parametrizes google executor parameters so that pipeline can now be run on different CloudOS environments (#281)
 - Adds a new `--download_from` option `FTP` mode to download SRA samples from [EBI FTP](https://ftp.sra.ebi.ac.uk/vol1/fastq/) (#283)
- Adds new parameter `--save_unmapped` that makes saving of STAR unmapped files optional (false by default) (#284)

#### Fixes:
 - Adds missing trimmomatic logs to the multiqc report (#244)
 - Implemented correct support for input strandness in star process when `--stranded` is `second-strand` (was hardcoded to `strType=2` and only supported `first-strand` or `false` before) (#264)
 - Issue that stringti_merged results folder as well as some other folders are missing all or some files (#263)
 - Fix pipeline crash when `params.stranded` was set to `false` (#276)
 - Fixes old parameters in google.config that were undesirably overwriting nextflow.config parameters on CloudOS (#281, [217e202](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/pull/245/commits/217e202cab3264c9d2d4cafe80b2476a2d837a85))
 
#### Updates:
 - Updates the following tools: (#248)
   - **STAR** `2.7.3` -> `2.7.9a` NOTE: Requires a new index! (updated in test profile)
   - **Samtools** `1.10` -> `1.13`
   - **StringTie** `2.1.3b` -> `2.1.7`
   - **Gffread** `0.11.7` -> `0.12.7`
   - multiqc `1.8` -> `1.11`
   - deeptools `3.4.0` -> `3.5.1`
   - bioconductor-rtracklayer `1.46.0` -> `1.52.0`
   - gffcompare `0.11.2` -> `0.12.6`
   - bedtools `2.29.2` -> `2.30.0`
   - sra-tools `2.10.8` -> `2.11.0`
   - pigz `2.3.4` -> `2.6.0`
   - gdc-client `1.5.0` -> `1.6.1`
 - Moves all containers to https://hub.docker.com/u/anczukowlab

#### Maintenance:
 - Considerably reduces number of basic redundant CI tests by removing completely the `max_retries` matrix and `push` from `on: [push, pull_request]`
 - Adds CI test for sra-downloading pipeline pathway (only supported with docker profile for now) (#253)

 
## v 1.0 - Initial pipeline release
