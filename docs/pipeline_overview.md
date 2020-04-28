# Pipeline overview

You can read more about what the pipeline does [here](../README.md#introduction)

The pipeline contains the following files/folders in rough order of importance:

- [`main.nf`](../main.nf) the main Nextflow script which executes all of the different processes when run
- [`nextflow.config`](../nextflow.config) the main configuration file - it includes defaults for all parameters, containers & profiles/links to other configs
- [`main.pbs`](../main.pbs) bash script to submit SLURM job
- [`README.md`](../README.md) high-level documentation for the pipeline including links to all of the documentation
- [`docs`](../docs) more detailed documentation for the pipeline
    - [`pipeline_overview.md`](pipeline_overview.md) this overview of the pipeline structure
    - [`run_on_cloudos.md`](run_on_cloudos.md) documentation for how to run the pipeline on Lifebit's platform CloudOS
    - [`run_on_sumner.md`](run_on_sumner.md) documentation for how to run the pipeline on JAX's HPC Sumner
    - [`usage.md`](usage.md) documentation for running the pipeline and parameters
- [`conf`](../conf) configuration files
    - [`examples`](../conf/examples) example configurations specific to an analysis and for testing (which specify the input parameters)
    - [`executors`](../conf/executors) configuration files for different infrastructure (eg for running on Sumner & over Google Cloud)
- [`examples`](../examples) contains examples data files
    - [`analyses`](../examples/analyses) examples data files for specific analyses
    - [`testdata`](../examples/testdata) example data files for testing
- [`containers`](../containers) contains `Dockerfile`s (instructions) to build the Docker images
- [`original_scripts`](../original_scripts) contains the original Bash & Nextflow scripts
- [`.gitignore`](../.gitignore) contains a list of files for git to ignore (not to commit to the repository)


The file structure looks as follows:
```
├── README.md
├── bin
│   ├── LU_postprocessing.R
│   ├── correct_gene_names.R
│   └── prepDE.py
├── conf
│   ├── examples
│   │   ├── MYC_MCF10A_0h_vs_MYC_MCF10A_8h.config
│   │   ├── human_test.config
│   │   └── test.config
│   └── executors
│       ├── base.config
│       ├── google.config
│       └── sumner.config
├── containers
│   ├── rmats4
│   │   ├── Dockerfile
│   │   └── environment.yml
│   └── splicing-pipelines-nf
│       ├── Dockerfile
│       ├── environment.yml
│       └── tagXSstrandedData.awk
├── docs
│   ├── github.md
│   ├── pipeline_overview.md
│   ├── run_on_cloudos.md
│   ├── run_on_sumner.md
│   └── usage.md
├── examples
│   ├── analyses
│   │   ├── MCF10_MYCER.datafiles.csv
│   │   ├── MYC_MCF10A_0h_vs_MYC_MCF10A_8h
│   │   │   ├── reads.csv
│   │   │   ├── reads_google_cloud.csv
│   │   │   └── rmats_pairs.txt
│   │   └── PRJNA453538.SraRunTable.txt
│   ├── assets
│   │   ├── TruSeq3-PE.fa
│   │   ├── TruSeq3-SE.fa
│   │   └── multiqc_config.yaml
│   └── testdata
│       ├── human_test
│       │   ├── get_human_paired_end_replicates.sh
│       │   ├── human_test_reps.csv
│       │   ├── human_test_reps_google.csv
│       │   └── rmats_pairs.txt
│       └── single_end
│           ├── rmats_pairs.txt
│           ├── test_reps.csv
│           └── test_reps_google.csv
├── main.nf
├── main.pbs
├── nextflow.config
└── original_scripts
    ├── bash
    │   ├── pipeline_splicing_with_arguments_parallel_part1.pbs
    │   ├── pipeline_splicing_with_arguments_parallel_part2.pbs
    │   ├── postprocessing.pbs
    │   └── run_pipeline_Olga_job1.sh
    ├── nextflow
    │   ├── main.nf
    │   └── rMATS_pipeline_samtools.nf
    └── postprocessing
        ├── LU_postprocessing.R
        ├── README.md
        ├── config_for_rmats_and_postprocessing.txt
        ├── config_for_splicing_pipeline_part2.txt
        ├── correct_gene_names.R
        ├── pipeline_splicing_with_arguments_parallel_part2_withGFFcompare.pbs
        ├── prepDE.py
        └── rmats_and_postprocessing.pbs
```