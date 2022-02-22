# Moving files from HPC to Cloud (particular to JAX Sumner)

#### Add singularity to $PATH:
module load singularity

#### Make some convenient commands to reduce typing:
alias gcloud="singularity exec /projects/researchit/crf/containers/gcp_sdk.sif gcloud"
alias gsutil="singularity exec /projects/researchit/crf/containers/gcp_sdk.sif gsutil"

#### Login to gcloud; this will return a url that you need to paste into a browser, which will take you through the google authentication process; you can use your jax email as userid and jax password to get in. Once you authenticate, it will display a code that you need to paste into the prompt provided in your ssh session on Sumner:
gcloud auth login --no-launch-browser

#### See which projects you have access to:
gcloud projects list

#### What is the project you are currently associated with:
gcloud config list project

#### Change project association:
gcloud config set project my-project

#### See what buckets are associated with my-project:
gsutil ls

#### See contents of a particular bucket:
gsutil ls -l gs://my-bucket

#### Recursively copy large directory from file system accessible on Sumner to your bucket:
gsutil -m -o GSUtil:parallel_composite_upload_threshold=150M cp -r my_dir gs://my_bucket/my_dir

#### Recursively copy a directory from your bucket to an existing directory on Sumner:
gsutil -m -o GSUtil:parallel_composite_upload_threshold=150M cp -r gs://my_bucket/my_dir my_dir
