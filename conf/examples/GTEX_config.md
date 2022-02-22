A minimal set of params need to run when downloading option is GTEX. Test is done with following params on a dev environment.

```yaml
params {
    reads = splicing-pipelines-nf/examples/GTEX/reads.csv
    manifest = manifest.json
    run_name = gtex_gen3
    download_from = GTEX
    key_file = credentials.json
    gtf = gencode.v32.primary_assembly.annotation.gtf
    star_index = /mnt/shared/gcp-user/session_data/star_75
    assembly_name = GRCh38
    readlength = 75
    stranded = false
    gc_disk_size = 200.GB
}
```
