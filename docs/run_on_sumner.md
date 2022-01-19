# Running your own analysis on Sumner

**Before running anything, make sure the pipeline is up to date by doing the following:**

Go into the splicing pipeline folder
```
cd /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/
```

Update local version of the pipeline. Note: you will need to enter your github username and password.
```
git pull
```

Note: if you have not successfully completed the pipeline test, see [here](../README.md##quick-start-on-sumner-jaxs-hpc)

This pipeline can be run on Sumner in three ways:
  1. Input a `reads.csv` file to input fastq files and run the pipeline in its entirety.
  2. Input a `reads.csv` file to input fastq files and run pipeline until STAR mapping step with `--test` parameter set to `true`.
  3. Input a `bams.csv` file to input bams and run steps of pipeline following STAR mapping (Stringtie and rMATS).

The `cacheDir` stores singularity images. This is set in splicing-pipelines-nf/conf/executors/sumner.config. For non-Anczukow users, this should be changed to a home directory.

## Running full pipeline with FASTQ input
### 1. Create a new run directory

All analyses should be run in `/projects/anczukow-lab/NGS_analysis`. If your dataset already has a folder here, use that directory. Otherwise, create a new directory with the same name found in the `/projects/anczukow-lab/fastq_files/` folder (Example DATASET Directory - `/projects/anczukow-lab/NGS_analysis/Dataset_4_MYC_MCF10A`). 

Create a new run directory within the appropriate dataset directory with the following format: runNumber_initials_date `run1_LU_20200519` (Example RUN Directory - `/projects/anczukow-lab/NGS_analysis/Dataset_4_MYC_MCF10A/run1_LU_20200519`).  

### 2. Create/Locate `reads.csv` file for your dataset

Input reads are specified by the `reads` input parameter, specifying a path to a CSV file. The format of CSV file will vary slightly based upon the data, see examples for:

- [single-end](../examples/testdata/single_end/test_reps.csv) - must contain columns for `sample_id` and `fastq`
- [paired-end](../examples/human_test/human_test_reps.csv) - must contain columns for `sample_id`, `fastq1` and `fastq2`

The 'reads.csv' column names must match the above [single-end] and [paired-end] examples. The `sample_id` can be anything, however each must be unique. The `fastq` column(s) should contain the path to FASTQ files (publicly accessible ftp, s3 and gs links are also accepted). You can create this on your local computer in excel and use WinSCP to move it to Sumner, or use create it using `nano` on the cluster.

There should be one `reads.csv` file per dataset. If your dataset already has a `reads.csv` file, proceed to step 2.


### 3. Create `rmats_pairs.txt` input file

Each rMATS comparison must be specified with a comparison name as well as the `sample_id` as specified in the [`reads`](../examples/testdata/human_test/human_test_reps.csv) file. See example [`rmats_pairs.txt`](../examples/human_test/rmats_pairs.txt). Each line in the file corresponds to an rMATS execution. The first column corresponds to a unique name/id for the rMATS comparison (this will be used for the output folder/file names)

* Replicates should be comma separated and the samples for the `b1` / `b2` files i.e. case and control should be space separated. b1 - control and b2 - case.
    <details>
    <summary>See examples</summary>

    #### Single sample pair:
    ```
    comparison_id[space]sample1[space]sample2
    ```

    #### Multiple sample pairs, no replicates:
    ```
    comparison1_id[space]sample1[space]sample2
    comparison2_id[space]sample3[space]sample4
    ```

    #### Multiple sample pairs, with multiple replicates:
    ```
    comparison1_id[space]sample1replicate1,sample1replicate2,sample1replicate3[space]sample2replicate1,sample2replicate2,sample2replicate3
    comparison2_id[space]sample3replicate1,sample3replicate2,sample3replicate3[space]sample4replicate1,sample4replicate1,sample4replicate1
    ```
    
     #### B1 only, no rMATS comparison (if this is run, set '--statoff' parameter to 'true'):
    ```
    comparison_id[space]sample1,sample2,sample3
    ```
    </details>


### 4. Setup `NF_splicing_pipeline.config`

This config file will be specific to your user and analysis. **You do not need to edit the pipeline code to configure the pipeline**. Descriptions of all possible parameters and their default values can be found [here](usage.md#all-available-parameters) and [here](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/03c977a4a5b386a1ea31c8aae78592432e38f3a2/nextflow.config). 

To create your own custom config (to specify your input parameters) you can copy and edit this [example config](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/d5ef6327fd3773484ecd1f0cdd593531e47cff39/NF_splicing_pipeline.config) file. 

**VERY IMPORTANT NOTES***

- Each time you run the pipeline, go through all possible parameters to ensure you are creating a config ideal for your data. If you do not specify a value for a parameter, the default will be used. All parameters used can be found in the `log` file. WHEN IN DOUBT, SPECIFY ALL PARAMETERS!

- You must name your config file `NF_splicing_pipeline.config` (as specified in main.pbs)

- Your `NF_splicing_pipeline.config` must be in the directory that you are running your analysis.

- The `readlength` here should be the length of the reads - if read length is not a multiple of 5 (ex- 76 or 151), set 'readlength' to nearest multiple of 5 (ex- 75 or 150). This extra base is an artifact of Illumina sequencing

- To run full pipeline, you **must** specify the following: `reads.csv`, `rmats_pairs.txt`, `readlength`, `assembly_name`, `star_index`, and `reference gtf`. This string can be a relative path from the directory in which you run Nextflow in, an absolute path or a link. 

- The star indexes must be generated prior to executing the pipeline (this is a separate step). 

- Currently, the two options for genomes are hg38 and mm10. If you wish to use a newer version of the genome, you will need to add this to the post-processing script.

### 5. Run the pipeline!

Ensure you have `NF_splicing_pipeline.config` in this directory. 

Run the pipeline! 
```
sbatch /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/main.pbs
```

## Running Stringtie and rMATS with BAM input
### 1. Create a new run directory

Create a new run directory within the appropriate dataset directory with the following format: runNumber_initials_date `run1_LU_20200519` (Example RUN Directory - `/projects/anczukow-lab/NGS_analysis/Dataset_4_MYC_MCF10A/run1_LU_20200519`).  

### 2. Create/Locate `bams.csv` file for your dataset

Input reads are specified by the `bams` input parameter, specifying a path to a CSV file.

- (create example) must contain columns for `sample_id`, `bam`, and `bam.bai`

The 'bams.csv' column names must match the above example. The `sample_id` can be anything, however each must be unique. The `bam` column should contain the path to BAM files. The `bam.bai` column should contain the path to BAM.BAI files. You can create this on your local computer in excel and use WinSCP to move it to Sumner, or use create it using `nano` on the cluster.

Supplying the `bams.csv` will signal to the pipeline to skip the first steps of the pipeline and start with Stringtie. No other parameter is needed.

### 3. Create `rmats_pairs.txt` input file

Each rMATS comparison must be specified with a comparison name as well as the `sample_id` as specified in the [`bams.csv`](create example) file. See example [`rmats_pairs.txt`](../examples/human_test/rmats_pairs.txt). Each line in the file corresponds to an rMATS execution. The first column corresponds to a unique name/id for the rMATS comparison (this will be used for the output folder/file names).

* Replicates should be comma separated and the samples for the `b1` / `b2` files i.e. case and control should be space separated
    <details>
    <summary>See examples</summary>

    #### Single sample pair:
    ```
    comparison_id[space]sample1[space]sample2
    ```

    #### Multiple sample pairs, no replicates:
    ```
    comparison1_id[space]sample1[space]sample2
    comparison2_id[space]sample3[space]sample4
    ```

    #### Multiple sample pairs, with multiple replicates:
    ```
    comparison1_id[space]sample1replicate1,sample1replicate2,sample1replicate3[space]sample2replicate1,sample2replicate2,sample2replicate3
    comparison2_id[space]sample3replicate1,sample3replicate2,sample3replicate3[space]sample4replicate1,sample4replicate1,sample4replicate1
    ```
    
     #### B1 only, no rMATS comparison (if this is run, set '--Statoff' parameter to 'true'):
    ```
    comparison_id[space]sample1,sample2,sample3
    ```
    </details>


### 4. Setup `NF_splicing_pipeline.config`

This config file will be specific to your user and analysis. **You do not need to edit the pipeline code to configure the pipeline**. Descriptions of all possible parameters and their default values can be found [here](usage.md#all-available-parameters) and [here](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/03c977a4a5b386a1ea31c8aae78592432e38f3a2/nextflow.config). 

To create your own custom config (to specify your input parameters) you can copy and edit this [example config]() file.  

**VERY IMPORTANT NOTES***

- Each time you run the pipeline, go through all possible parameters to ensure you are creating a config ideal for your data. If you do not specify a value for a parameter, the default will be used. All parameters used can be found in the `log` file. WHEN IN DOUBT, SPECIFY ALL PARAMETERS!

- You must name your config file `NF_splicing_pipeline.config` (as specified in main.pbs).

- Your `NF_splicing_pipeline.config` must be in the directory that you are running your analysis.

- The `readlength` here should be the length of the reads - if read length is not a multiple of 5 (ex- 76 or 151), set 'readlength' to nearest multiple of 5 (ex- 75 or 150). This extra base is an artifact of Illumina sequencing

- Currently, the two options for genomes are hg38 and mm10. If you wish to use a newer version of the genome, you will need to add this to the post-processing script.

### 5. Run the pipeline!

Ensure you have `NF_splicing_pipeline.config` in this directory. 

Run the pipeline! 
```
sbatch /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/main.pbs
```
