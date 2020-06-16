# Containers introduction

Containers can be used to bundle software dependencies. This means that analyses can be made more reproducible, installation of tools much easier and therefore that the same workflow or analysis can be run across multiple different compute infrastructures in a very portable manner. For example, this workflow can be run on HPC (using [JAX](https://www.jax.org/)'s Sumner HPC) or on the cloud over ([Lifebit's CloudOS](https://lifebit.ai/cloudos) platform with AWS & GCloud). This is thanks to containers and also due to the workflow manager [Nextflow](https://www.nextflow.io) which has in-built support containers such as Docker and Singularity.

## Resources
- An introduction to containers (as well as Nextflow & CloudOS) can be found [here](https://github.com/lifebit-ai/jax-tutorial#session-2-docker)
- Instructions on installing Docker and Singularity can be found [here](https://github.com/lifebit-ai/jax-tutorial/blob/master/README.md#ii-installing-docker) :warning: You will need root permissions to install Docker :warning:
- You can also read this [this guide](https://docs.docker.com/get-started/) (4mins read) for a more high level overview of Docker containers.

If your still confused after reading it don't hesitate to ping [@PhilPalmer](https://github.com/PhilPalmer) or [@adeslatt](https://github.com/adeslatt). 

## Important note

One important note is that containers would ideally be:
1. Docker containers
    - This is because Nextflow can convert Docker -> Singulairty containers but not the other way around
2. Hosted in a [google container registry (gcr)](https://cloud.google.com/container-registry)
    - This is because to run on Google Cloud (eg on CloudOS) containers will take too long to be fetch and cause the pipeline to fail if they are not hosted here

Doing both of these things makes containers as portable as possible. However, both also have issues because Docker requires root or admin acess (which is not available on Sumner) & you must have access to a gcr to be able to push containers there.

## To build a Docker container
If you need to modify one of the containers, eg to update or add more software dependencies you can do so like so:
```bash
docker build -t <registry_user>/<image_name>:<tag> .
```

Eg:
```
cd containers/splicing-pipelines-nf/
docker build -t gcr.io/nextflow-250616/splicing-pipelines-nf:gawk .
```

## To push your image to a container register
```bash
docker push <registry_user>/<image_name>:<tag>
```

### Eg to GCR
:warning: You will need credentials/access for a google container registry :warning:
```
gcloud auth login
docker push gcr.io/<project_id>/<image_name>:<tag>
```
