#!/bin/bash
samplesheet_path="/home/ubuntu/data/local/MD_project/scripts/pipelines/exome/04_samplesheets/samplesheets/all_2.csv"
JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/21.0.2-tem"
JAVA_CMD="/home/ubuntu/.sdkman/candidates/java/21.0.2-tem"
java -version


nextflow run nf-core/sarek \
--input $samplesheet_path \
-profile singularity \
-c /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/05_sarek/sarek.config \
--fasta "/data/local/MD_project/data/reference/human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna" --fasta_fai "/home/ubuntu/data/local/MD_project/data/reference/human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.fai" \
--dict "/home/ubuntu/data/local/MD_project/data/reference/human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.dict" \
--tools mutect2,haplotypecaller,ascat,controlfreec,msisensorpro --wes \
--outdir "/home/ubuntu/data/local/MD_project/data/exome/processed/04_sarek"  -w "/home/ubuntu/data/local/MD_project/scripts/pipelines/exome/05_sarek/work2" -resume

