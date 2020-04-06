# Run on CloudOS

## 0) Create an account & join a team
[Sign-up to CloudOS](https://deploit.lifebit.ai/register) if you haven't already and join your team workspace. (While you don't have to join a team it is important to have your cloud account linked. Joining your team means that you won't need to do this and can collaborate with colleagues).

## 1) Import the pipeline

If you navigate to the pipelines page (this can be found in the navigation bar on the left-hand-side) and then can see the pipeline under `MY PIPELINES & TOOLS` continue to step 2

Otherwise, you will need to import it.

To import the pipeline you need to:
- Go to the pipelines page
- Click the green `New` button (top right)
- `Select` the GitHub icon to import the Nextflow pipeline from GitHub
- Paste the URL of our pipeline [`https://github.com/lifebit-ai/splicing-pipelines-nf`](https://github.com/lifebit-ai/splicing-pipelines-nf)
- Name the pipeline, eg `splicing-pipeline`
- (Optional:) enter a pipeline description
- Click `Next` & `Create pipeline` :tada:

![import_splicing_pip](https://raw.githubusercontent.com/lifebit-ai/images/master/jax_splicing/import_splicing_pip.gif)

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
nextflow run https://github.com/PhilPalmer/splicing-pipelines-nf
--reads 'https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/pipelines/rmats/reads/replicates/paired_end_reads_replicates.csv'
--genome GRCh38
--b1 'https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/pipelines/rmats/reads/replicates/b1.txt'
--b2 'https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/pipelines/rmats/reads/replicates/b2.txt'
--overhang 100
--igenomes_base 's3://ngi-igenomes/igenomes/'
--gencode_gtf 'https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/jax/splicing-pipelines-nf/references/GRCh38/gencode.v32.annotation.gtf'
--max_cpus 16
--max_memory '120.GB'
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

This page is completely sharable, for example, you can view a successfully completed example job [here](https://cloudos.lifebit.ai/public/jobs/5e831015e7d1990104cb8090) 
![splicing_pip_job_page](https://raw.githubusercontent.com/lifebit-ai/images/master/jax_splicing/splicing_pip_job_page.png)