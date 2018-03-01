#!/bin/bash

# Load modules
module load bioinfo-tools \
    sratools/2.8.0 \
    Salmon/0.8.2 \
    star/2.5.3a \
    samtools/1.5 \
    bamtools/2.3.0 \
    picard/2.10.3 \
    GATK/3.8-0 \
    snpEff/4.2

# Loop over single and paired-end reads
for CURRENT_LAYOUT in SINGLE PAIRED; do

    # Run snakemake
    echo "Processing series with ${CURRENT_LAYOUT}-END reads"
    snakemake \
        --snakefile Snakefile \
        --config LAYOUT=$CURRENT_LAYOUT \
        --config group_col="group" \
        --jobs 20 \
        --keep-going \
        --cluster-config cluster.yaml \
        --cluster "sbatch \
                      --account {cluster.account} \
                      --time {cluster.time} \
                      --partition {cluster.partition} \
                      --ntasks {cluster.ntasks} \
                      --mail-type {cluster.mail_type} \
                      --job-name {cluster.job_name} \
                      --output=/dev/null \
                      --error=/dev/null"
done

# Send notification mail
echo "" | mail -s "Snakemake finished." erikfas@kth.se
