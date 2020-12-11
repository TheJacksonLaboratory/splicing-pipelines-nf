A minimal set of params need to run when downloading option is GEN3-DRS. Test is done with following params on a dev environment.

```yaml
params {
    reads = splicing-pipelines-nf/examples/gen3/manifest.csv
    run_name = gtex_gen3
    download_from = GEN3-DRS
    key_file = gen3_cred.json
    gtf = gencode.v32.primary_assembly.annotation.gtf
    star_index = /mnt/shared/gcp-user/session_data/star_80
    genome_fasta = s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/genome.fa
    assembly_name = GRCh38
    readlength = 80
    statoff = true
    max_cpus = 14
    max_memory = 60
}
```
