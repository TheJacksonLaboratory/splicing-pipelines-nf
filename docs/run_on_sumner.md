# Test the pipeline

See [here](../README.md##quick-start-on-sumner-jaxs-hpc)

# Running your own analysis

## 1. Create `reads` file 

* This file contains a `sample_id` and path to fastq files for each sample 

* You will need to create a CSV file containing the path to your input `reads`. You can see examples for [single-end](../examples/testdata/single_end/test_reps.csv) and [paired-end](../examples/human_test/human_test_reps.csv) data

* These files must have the column names as in the above examples

* The `sample_id` can be anything, however each must be unique

* You can create this on your local computer in excel and use WinSCP to move it to Sumner, or use create it using `nano` on the cluster.

### If you wish to run rMATS you will need to create `rmats_pairs` input file

* Each rMATS comparison must be specified with a comparison name as well as the `sample_id` as specified in the [`reads`](../examples/testdata/human_test/human_test_reps.csv) file. See example [`rmats_pairs.txt`](../examples/human_test/rmats_pairs.txt)

* Each line in the file corresponds to a rMATS execution

* The first column corresponds to a unique name/id for the rMATS comparison (this will be used for the output folder/file names)

* Replicates should be comma seperated and the samples for the `b1` / `b2` files i.e. case and control should be space seperated
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
    </details>


## 2. Setup your own configuration file

Here you will be [adding your own custom config](https://nf-co.re/usage/configuration#custom-configuration-files). The config file will be specific to your user and analysis. **You do not need to edit the pipeline code to configure the pipeline**.

* To create your own custom config (to specify your input parameters) you can copy and edit this [example config](../conf/examples/MYC_MCF10A_0h_vs_MYC_MCF10A_8h.config) file.

* **You must name your config file `NF_splicing_pipeline.config`**

* Find more information on all available parameters [here](usage.md#all-available-parameters) and see their default values in [`nextflow.config`](../nextflow.config). **NOTE**: You do not need to specify all parameters if the default parameter is acceptable

* You will need to specify the path to your `reads` and `rmats_pairs` input files. This string can be a relative path from the directory which you run Nextflow in, an absolute path or even a link.

## 3. Run the pipeline

* If you have not done so already, create a new run directory within the appropriate dataset directory with the following format: runNumber_initials_date `run1_LU_20200519`

* Ensure you have the following files in this directory: `reads.csv`, `rmats_pairs.txt`, `NF_splicing_pipeline.config`

* Run the pipeline! 
```
sbatch /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/main.pbs
```
**NOTE: if running on cloud we don't have main.pbs set up yet**


# Below has not been updated yet (5/19/20)
## Examples

See [`MYC_MCF10A_0h_vs_MYC_MCF10A_8h.config`](../conf/examples/MYC_MCF10A_0h_vs_MYC_MCF10A_8h.config) for an example analysis comparing the 0h and 8h timepoints

The analysis compares `MYC_MCF10A_0h` with 3 replicates and `MYC_MCF10_8h`.
The details of what needs to be configured to do this comparison analysis are found in three files:

All the analyses can be kept in [`examples/analyses`](../examples/analyses) subdirectory (by convention). Encoded in the file name is the metadata details that outline the comparison that is being completed.  In this case capturing the statement above (`MYC_MCF10A_0h` vs `MYC_MCF10A_8h`).

These files can be specified via command line or via a config file.

To specify via command line:

* `--reads examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/reads.csv`
    This file contains the `sample_id` a short name uniquely defines the sample within this comparison
    comma seperated with the complete path for the `left` and `right` `fastqs`.   
    
* `--rmats_pairs examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/rmats_pairs.txt`
    This is a space/comma separated file containing the a list of the rMATS comparisons you wish to perform.

# Bonus: useful Nextflow options

Whereas parameters are set on the command-line using double dash options eg `--reads`, parameters passed to Nextflow itself can be provided with single-dash options eg `-profile`.

You can see some of these options [here](https://www.nextflow.io/docs/latest/tracing.html) in the Nextflow documentation.

Some useful ones include:
- `-resume` which will [resume](https://www.nextflow.io/docs/latest/getstarted.html?highlight=resume#modify-and-resume) any cached processes that have not been changed
- `-with-trace` eg `-with-trace trace.txt` which gives a [trace report](https://www.nextflow.io/docs/latest/tracing.html?highlight=dag#trace-report) for resource consumption by the pipeline
- `-with-dag` eg `-with-dag flowchart.png` which produces the [DAG visualisation](https://www.nextflow.io/docs/latest/tracing.html?highlight=dag#dag-visualisation) graph showing each of the different processes and the connections between them (the channels)
