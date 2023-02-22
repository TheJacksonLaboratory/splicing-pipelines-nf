# Running your own analysis on GCP

**Make sure you have the gcp-splice branch, it is the only branch tested for running the pipeline directly on GCP:**

Go into the splicing pipeline folder
```
cd /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/
git branch -a #### to list all brances
git checkout gcp-splice #### to checkout in branch testing for running the pipeline on GCP
git branch #### to check which branch you are on
```

The steps to run the nextflow splicing pipeline are listed below:

## Running full pipeline with FASTQ input

### 1. Create conda environment
```
conda env list
conda env create --name splicing-pipelines-nf -f containers/splicing-pipelines-nf/environment.yml
conda activate splicing-pipelines-nf
```
### 2. Install nextflow in the conda env
```
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
```

### 3. Install google SDK on your local and make sure you have access to jax-cloudos-olgaanczukow project
It is good to follow the instruction provided at (https://cloud.google.com/sdk/docs/install) .... Once the installation is done
```
gcloud init
gcloud auth login
gcloud auth application-default login
gcloud config configurations list
gcloud config set project jax-cloudos-olgaanczukow
```

### 4. Create/Locate `reads.csv` file for your dataset

Input reads are specified by the `reads` input parameter, specifying a path to a CSV file. The format of CSV file will vary slightly based upon the data, see examples for:

- [paired-end](jax-cloudos-olgaanczukow-project-data/test_data/gcp-test-humanReps.csv) - must contain columns for `sample_id`, `fastq1` and `fastq2`

### 5. Create `rmats_pairs.txt` input file

Each rMATS comparison must be specified with a comparison name as well as the `sample_id` as specified in the [`reads`](jax-cloudos-olgaanczukow-project-data/test_data/gcp-test-humanReps.csv) file. See example [`rmats_pairs.txt`]([../examples/human_test/rmats_pairs.txt](https://storage.cloud.google.com/jax-cloudos-olgaanczukow-project-data/test_data/rmats_pairs.txt). Each line in the file corresponds to an rMATS execution. The first column corresponds to a unique name/id for the rMATS comparison (this will be used for the output folder/file names)

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


### 4. Setup using the correct configuration file

- To run the analysis specifically on GCP the most important configuration file is: `conf/executors/google.config`

**You might have to edit the configurtion file occassionally to specify a particular VM**. Descriptions of all possible parameters and their default values can be found [here](usage.md#all-available-parameters) and [here](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/03c977a4a5b386a1ea31c8aae78592432e38f3a2/nextflow.config). 

- Another important configuration file `nextflow.config` 
Make sure the flag allowing nextflow to run on GCP is present in the file and the flags for GCP are passed correctly in the configuration file
```
google { includeConfig 'conf/executors/google.config' }
```
**VERY IMPORTANT NOTES***

- To run full pipeline, you **must** specify the following: `reads.csv`, `rmats_pairs.txt`, `readlength`, `assembly_name`, `star_index`, and `reference gtf`. This string can be a relative path from the directory in which you run Nextflow in, an absolute path or a link. 

### 5. Plcing the dataset on GCP bucket
It is more convenient to place all the datasets, including `fastq files, read.csv, rmats.txt, annotation.gtf, reference fasta and star alignment` in the bucket on GCP. For example please look at:
  - (jax-cloudos-olgaanczukow-project-data/human_reference_genome) .... for storing files related to reference
  - (jax-cloudos-olgaanczukow-project-data/test_data) ... for storing the dataset files

### 6. Run the pipeline!

Ensure you have `nextflow.config and NF_splicing_pipeline.config` in running directory. 

Run the pipeline! 
```
nextflow run main.nf -c nextflow.config -work-dir gs://jax-cloudos-olgaanczukow-project-data/test_data/splicing-work-dir/output -profile google -resume
```
