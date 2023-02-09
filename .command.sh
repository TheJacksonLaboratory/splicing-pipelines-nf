#!/bin/bash -ue
fastqc --casava --threads 4 SRR12924637_1.fastq.gz SRR12924637_2.fastq.gz

# save .command.* logs
task_hash=`basename ${PWD} | cut -c1-6`; mkdir command-logs-$task_hash ; cp .command.*{err,log,sh} command-logs-$task_hash
