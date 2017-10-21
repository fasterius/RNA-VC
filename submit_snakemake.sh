#!/bin/bash

# Run snakemake
snakemake \
    -s snakefile \
    --keep-going \
    --cluster-config cluster.yaml \
    --cluster "sbatch \
                  --account {cluster.account} \
                  --time {cluster.time} \
                  --partition {cluster.partition} \
                  --n-tasks {cluster.n_tasks} \
                  --job-name {cluster.job_name} \
                  --mail-type {cluster.mail_type}"