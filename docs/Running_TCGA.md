# This is how Laura ran TCGA 

## 1) Get file UUISs for tumors

You will need the file UUID for each tumor (ex: 2bf57145-58f2-4cbb-aa18-176caa1bd4fe)

For each tumor type, generate a .csv file containing all file IDs and read lengths. Inside `TCGA_and_GTEX_analysis_LU` folder, there is a subfolder for each tumor type. There you will find a script to get UUIDs (ex: `get_LAML_UUIDs.R`). The output from this script has already been saved for each tumor type in the `TCGA` folder (ex: `TCGA-LAML_UUIDs.xlsx`).

Note: this file contains extra information, but what is important is the file ID and read length. **DO NOT TRUST THE PAIRED END DOCUMENTATION. It is incorrect** 

## 2) Generate a reads.txt 

Use the UUID excel file to obtain the samples you wish to run. For tumor types with more than 100 samples, run each in a batch of 100. Generate a reads.txt file for each batch. [example](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/master/examples/TCGA/reads.txt)

## 3) Get GDC key file 

Naviate to the [GDC portal](https://portal.gdc.cancer.gov/) and log in. 
Click on your user name and download token. 

## 4) Upload reads.txt files to CloudOS

Refer [here](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/master/docs/run_on_cloudos.md) for infomration on how to upload. 

## 5) Run mapping step for TCGA 

If an old job exists, you can clone that job. However, if you wish to set up a new analysis, follow directions [here](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/master/docs/run_on_cloudos.md) but keep in mind the few parameters below: 

```
download_from tcga
sliding_window false 
key_file 'path to key file'
readlength 50
mismatch 5
filterScore 0.33
sjdbOverhangMin 1
test true
gc_disk_size 200.GB
```
Note: to run only mapping make sure to specify `test` as true. Also, depending on the number of samples you may wish to increase/decrease gc_disk_size 

Make sure to increase added storage when you select the instance. For mapping 100 samples, typically use 2000. 

## 6) Monitor each job, document run info and failed samples

Once job is completed, make sure to document run name, job ID, time, cost, and results path. 

Identify failed samples by clicking on the number below `failed` in the Status box. The failed sample and the process it failed on are indicated on the left. 

## 7) Re-run failed samples 

Samples can fail for a variety of reasons. Often, when you rerun them, they will pass. However, if they fail a second time, they will likely always fail. 

Generate a new reads.txt file containing file UUIDs for failed samples. Run these like you did above. Note which ones fail a second time. 

## 8) Get paths to bams and generate bams.csv for rMATS

- Obtain path to results folder for each batch 
- Write out path to bams in a text file `gsutil ls results_path/star_mapped/*/*.bam > batch1_bams.txt` This will give you a text file with the paths to all bams for that batch. 
- Open the bams.txt file in excel
- Make sure the number of paths is equal to the number of passed samples in that batch. If the number is off, there may be an issue with the results directory where some bams have different paths. Sometimes if a sample failed on Star, a bam file will still be generated and this path should be removed from your list. 
- Combine all bam paths from all batches into one excel document. Note the number of paths should equal the total number of passed samples 
- Insert blank column to the left of the bam paths
- Copy the column of paths and paste in a new column (example column D). 
- Highlight the copied paths and select Data and text to columns. Make sure delineated is selected and click next. Check the `other` box and type `/` in the text box. Click next and then finish. This will separate the path by `/`. 
- Copy column after star_mapped. This column has the UUID for each sample. Paste this column in the first column of the spreadsheet. 
- Delete all the columns you just split 
- Generate path to bam.bai in column C using the SUBSTITUTE function `SUBSTITUTE(B1, ".bam", ".bam.bai")`. Copy this all the way down the document. 
- Add new row at the top and specify column names `sample_id`, `bam`, `bam_index`. 
- Copy and paste as text into new excel document (for example, into run_info document). 

Note: if more than 100 samples, you will need to run rMATS in batches. Assign samples to batches of approximately equal sizes. These may not be the same batches you used in the mapping step. 

## 9) Generate bams.csv and rmats_pairs.txt

- Use the file generated above to make your `bams.csv`. Make sure you have column named `sample_id`, `bam`, `bam_index`. Note, you need one bams.csv per rMATS run
- For each rMATS run you will need an rMATS_pairs.csv. 
- Upload `bams.csv` and `rmats_pairs.txt` to CloudOS into the appropriate folders (ex: gs://anczukow-bucket/TCGA-GBM)

## 10) Run rMATS

- If you have a job you wish to clone, you may do so. Otherwise you may create a new analysis. Keep in mind the following parameters

```
single_end
readlength
statoff
gc_disk_size 2000.GB
max_retries 10 
```

NOTE: if you are running b1 only, set  `statoff` to true. Again you may need to change the gc_disk_size if larger datasets. 

- use the same n1-standard-8 instance. 
- Make sure to add storage. For example, rMATS for ~100 samples, 10000. For fewer samples you can decrease this number

## 11) Download rMATS files 

Download files you need as you would for any other dataset. 



