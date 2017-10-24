#!/bin/bash

# Loop over single and paired-end reads
for CURRENT_LAYOUT in SINGLE PAIRED; do

    # Run snakemake
    snakemake \
        --snakefile snakefile \
        --config LAYOUT=$CURRENT_LAYOUT \
        --jobs 20 \
        --keep-going \
        --cluster-config cluster.yaml \
        --cluster "sbatch \
                      --account {cluster.account} \
                      --time {cluster.time} \
                      --partition {cluster.partition} \
                      --ntasks {cluster.ntasks} \
                      --mail-type {cluster.mail_type} \
                      --job-name {cluster.job_name}"
done
