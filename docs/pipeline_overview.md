# Pipeline overview

You can read more about what the pipeline does [here](../README.md#introduction)

The pipeline contains the following files/folders of interest:

- `main.nf` the main Nextflow script
- `nextflow.config` the main configuration file - it includes defaults for all parameters, containers & links to other configs/profiles
- `main.pbs` bash script to submit SLURM job
- `README.md` high-level documentation for the pipeline with links
- `docs` more detailed documentation for the pipeline
- `examples` contains example testdata & analyses

There are also the following folders which can (for the most part be ignored):
- `conf` contains the configuration file
  * `base.config` to specify resources
  * `igenomes.config` specifying the path to the igenomes resources available
      To mirror the setup on the CloudOS Universal Research Platform, adopting this standard setup.
      * Currently in `/projects/adeslat/igenomes/` `GRCh38` has been downloaded
      * `STARIndex` downloaded from the `AWS S3 Budget` Amazon has donated to host these standards
- `containers` contains `Dockerfile`s/instructions to build the Docker images
- `original_scripts` contains the original Bash & Nextflow scripts
- `.gitignore` contains a list of files for git to ignore (not to commit to the repository)