#!/bin/bash
samplesheet_path="/home/ubuntu/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/bbsplit_test.csv"
java --version

nextflow run nf-core/sarek --input $samplesheet_path \
 -profile singularity --fasta "/home/ubuntu/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
--tools mutect2,haplotypecaller,ascat,controlfreec,msisensorpro --wes \
--outdir "/data/local/MD_scholarly/data/processed/exome/sarek" -c sarek.config


