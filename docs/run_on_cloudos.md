# Run on CloudOS

## 0) Create an account & join a team
[Sign-up to CloudOS](https://deploit.lifebit.ai/register) if you haven't already and join your team workspace. (While you don't have to join a team it is important to have your cloud account linked. Joining your team means that you won't need to do this and can collaborate with colleagues).

## 2) Run the pipeline

Pipelines can be run in three simple steps:
1. Select the pipeline
2. Select data & parameters
3. Run the analysis

### 2.1) Select the pipeline

Once the pipeline is imported it will automatically be selected.

Alternatively, you can navigate to the pipelines page. Where you can find the imported pipeline under `MY PIPELINES & TOOLS`. To select the pipeline you need to click the card for the pipeline.

### 2.2) Select data & parameters

Next we want to link to data & add parameters. You can navigate/import data on CloudOS or simply paste links to data

For example, for this pipeline you may want to add the following parameters and example data:
```
nextflow run https://github.com/lifebit-ai/splicing-pipelines-nf
--reads 'https://github.com/lifebit-ai/splicing-pipelines-nf/raw/master/examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/reads_google_cloud.csv'
--b1 'https://github.com/lifebit-ai/splicing-pipelines-nf/raw/master/examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/b1.txt'
--b2 'https://github.com/lifebit-ai/splicing-pipelines-nf/raw/master/examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/b2.txt'
--star_index 'gs://cloudosinputdata/inputs/splicing-pipelines-nf/references/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex'
--gtf 'gs://cloudosinputdata/inputs/splicing-pipelines-nf/gencode.v32.primary_assembly.annotation.gtf'
```

![run_splicing_pip](https://raw.githubusercontent.com/lifebit-ai/images/master/jax_splicing/run_splicing_pip.gif)

### 2.3) Run the analysis

Select a project & instance:

Before running the job you must:

1. Select (and possible create) the project (which is like a folder used to group multiple analyses/jobs), eg `Splicing`
2. Choose the instance to set the compute resources such as CPUs & memory. For example, `r4.4xlarge` on AWS
3. Finally, click Run job

## 3) Monitor the analysis
To monitor jobs you can click on the row for any given job. Immediately after running a job its status will be initialising. This normally occurs for ~5mins before you are able to view the progress of the job.

Once on the job monitor page, you can see the progress of the job update in real time. Information such as the resources i.e. memory & CPUs is displayed. Once the job has finished the results can be found in the results tab as well as any reports for select pipelines.

This page is completely sharable, for example, you can view a successfully completed example job [here](https://cloudos.lifebit.ai/public/jobs/5e87ef928079200103b0a0b8) 
![splicing_pip_job_page](https://raw.githubusercontent.com/lifebit-ai/images/master/jax_splicing/splicing_pip_job_page.png)

### Important notes while running large number of samples on CloudOS

- For one sample to run the first part of the pipeline (having steps - get_tcga_bams, bamtofastq, fastqc, trmmomatic, fastqc_trimmed, star, multiqc) require about 20 GB of results space per sample. So in order to accommodate that you need to specify 20*Number of samples in the CloudOS GUI while selecting an instance.  (Example 100 samples 20*100 = 2000 GB).

- Also accordingly there is a special parameter `--gc_disk_size` which is specific to google-cloud executor. Adds disk-space for few aggregative processes. (default: "200 GB" based on 100 samples. Which is 2GB x Number of Samples)

###Helpful Tips
[Import the Pipeline](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/master/docs/import_pipeline) 
