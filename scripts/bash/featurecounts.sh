#!/bin/bash

# Go to directory with BAM files
cd /mnt/c/Users/Nitin/Bulk_tutorial

# Make output quants folder
mkdir -p /mnt/c/Users/Nitin/Bulk_tutorial/quants

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 -a /mnt/c/Users/Nitin/Bulk_tutorial/Homo_sapiens.GRCh38.114.gtf \
        -o /mnt/c/Users/Nitin/Bulk_tutorialquants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "✅ Completed $bam in $runtime minutes."
    echo "------------------------------------"
done
