# Troubleshooting pipeline 

If you are getting errors when running any steps of the pipeline, here are some tips to help narrow down the issue 

## 1. Check `splicing.[job_number].err` file 

To view use `less` or `cat`, depending on file size. If there are errors listed in this file, try and address them. If this file is empty, move onto the next step. 

## 2. Check `splicing.[job_number].log` file 

This file contains a lot of useful information. At th ebeginning of the file you will see where it downloaded nextflow and set up the parameters for the pipeline. Check and make sure all of this information is correct. 

The main body of this file contains a list of steps in the pipeline and a progress report. It will indicate which step and which file are being processed. If the processes run successfully for all samples this will be indicated by X/X and a check mark. Scroll through this file and see if there are any errors. If there are, you can check the work directory for that step by using the path indicated on left of process name. 

## 3. Check the `nextflow.log`

This is a hidden file that can be accessed see using `ls -la` 
