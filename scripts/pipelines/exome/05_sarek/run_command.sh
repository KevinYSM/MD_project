#!/bin/bash
samplesheet_path="/home/ubuntu/data/local/MD_project/scripts/pipelines/exome/04_samplesheets/samplesheets/all.csv"
java --version
JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/current/bin/java"

nextflow run nf-core/sarek --input $samplesheet_path \
 -profile singularity --fasta "/data/local/MD_project/data/reference/human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna" \
--tools mutect2,haplotypecaller,ascat,controlfreec,msisensorpro --wes \
--outdir "/home/ubuntu/data/local/MD_project/data/exome/processed/04_sarek" -c /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/05_sarek/sarek.config


