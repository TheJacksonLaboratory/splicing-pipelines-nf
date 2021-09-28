# Changelog

## v 1.1 - Pipeline improvements

### Fixes:
 - Added missing trimmomatic logs to the multiqc report
 - Implemented correct support for input strandness in star process when `--stranded` is `second-strand` (was hardcoded to `strType=2` and only supported `first-strand` or `false` before)

### Updates:
 - Updated the following tools:
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
 - Moved all containers to https://hub.docker.com/u/anczukowlab

### Maintenance:
 - Consideably reduced number of basic redundant CI tests by removing completely the `max_retries` matrix and `push` from `on: [push, pull_request]`
 - Added CI test for sra-downloading pipeline pathway (only supported with docker profile for now)

### Enhancements:
 - Added saving of all the process .command* log files to results/process-logs folder
 - Added pipeline workdir `--cleanup` option to clear all intermediate files on pipeline successful completion
 - Added pipeline `--error_strategy` parameter to be able to specify pipeline error strategy directly from command line (doesn't work if specified in config linked by `-c` or `-config` nextflow params)

## v 1.0 - Initial pipeline release