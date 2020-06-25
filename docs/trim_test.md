# Run trimming and mapping test 

This test is used to determine the optimal read length if user wishes to crop all reads to one length. This test will run FastQC, trimmomatic, and Star on a subset of your data at a variety of read lengths. Things to consider when deciding: 

* Percentage of reads remaining after trimmomatic

* Star mapping stats

* Longer reads produce higher quality data for splicing analysis

If you have not done so already, create a new run directory within the appropriate dataset directory with the following format: runNumber_initials_date `run1_LU_20200519`


## 1. Create `NF_splicing_pipeline.config`

For example see [here](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/master/conf/examples/trim_test.config). **You must name your config file `NF_splicing_pipeline.config`**

For `readlength`, enter the length of the original reads (as multiple of 5)

* Example: 75, 80, 90, etc.

## 2. Choose an `increment`

What increments of read lengths would you like to test? For example, if my original read length is 150, I might want to test increments of 10. Meaning 140bp, 130bp, 120bp, and 110bp. 
    
Note: This test will only test 4 different read lengths. If you wish to test more, rerun the test and decrease the value specified for full length `readlength`

## 3. Run `trim_test.pbs`

To run: 
```
sbatch --export=ALL,increment=10 /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/trim_test.pbs
```
Note: Replace `increment` value with your desired value

## 4. Use MulitQC output to determine optimal read length

Copy the multiQC reports for each read length to your local computer using WinSCP. Use these to look at the trimming and mapping stats. Alternatively you can look at the log files generated. 

NOTE: all multiQC files have the same name. You will need to either change their name or put them in different folders on your computer. 
