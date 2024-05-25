#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
log_file=logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}')."$RANDOM".log
exec &> >(tee -a "$log_file")

echo "$CONDA_PREFIX"

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/21.0.2-tem"
JAVA_CMD="/home/ubuntu/.sdkman/candidates/java/21.0.2-tem"
SINGULARITY_TMPDIR="/data/local/MD_project/scripts/pipelines/exome/05_sarek_joakim/pipeline_tmp/s_tmp"
SINGULARITY_LOCALCACHEDIR="/data/local/MD_project/scripts/pipelines/exome/05_sarek_joakim/pipeline_tmp/s_cache"
#Works for Sarah:
#Nextflow: 22.04.3
#Singularity: 3.8.6
#Sarek: 3.1.2

nextflow run nf-core/sarek \
    -profile singularity \
    -c custom.config \
    -w work/batch_1 \
    -resume \
    --intervals /data/local/reference/agilent/SureSelect_XT_HS_Human_All_Exon_V8/hg38/S33266340_Padded.reformatted.bed \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --max_memory '120.GB' \
    --input inputs/batch_1.csv \
    --outdir /data/local/MD_project/data/exome/processed/04_sarek_joakim \
    --genome GATK.GRCh38 \
    --igenomes_base /data/local/reference/aws/igenomes \
    --max_cpus 63 \
    --tools mutect2,controlfreec,cnvkit,msisensorpro,vep,tiddit,haplotypecaller,manta,snpeff \
    --wes \
    --tmp-dir /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/05_sarek_joakim/pipeline_tmp/s_tmp


nextflow run nf-core/sarek \
    -profile singularity \
    -c custom.config \
    -w work/batch_2 \
    -resume \
    --intervals /data/local/reference/agilent/SureSelect_XT_HS_Human_All_Exon_V8/hg38/S33266340_Padded.reformatted.bed \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --max_memory '120.GB' \
    --input inputs/batch_2.csv \
    --outdir /data/local/MD_project/data/exome/processed/04_sarek_joakim \
    --genome GATK.GRCh38 \
    --igenomes_base /data/local/reference/aws/igenomes \
    --max_cpus 63 \
    --tools mutect2,controlfreec,cnvkit,msisensorpro,vep,tiddit,haplotypecaller,manta,snpeff \
    --wes


nextflow run nf-core/sarek \
    -profile singularity \
    -c custom.config \
    -w work/batch_3 \
    -resume \
    --intervals /data/local/reference/agilent/SureSelect_XT_HS_Human_All_Exon_V8/hg38/S33266340_Padded.reformatted.bed \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --max_memory '120.GB' \
    --input inputs/batch_3.csv \
    --outdir /data/local/MD_project/data/exome/processed/04_sarek_joakim \
    --genome GATK.GRCh38 \
    --igenomes_base /data/local/reference/aws/igenomes \
    --max_cpus 63 \
    --tools mutect2,controlfreec,cnvkit,msisensorpro,vep,tiddit,haplotypecaller,manta,snpeff \
    --wes


nextflow run nf-core/sarek \
    -profile singularity \
    -c custom.config \
    -w work/batch_4 \
    -resume \
    --intervals /data/local/reference/agilent/SureSelect_XT_HS_Human_All_Exon_V8/hg38/S33266340_Padded.reformatted.bed \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --max_memory '120.GB' \
    --input inputs/batch_4.csv \
    --outdir /data/local/MD_project/data/exome/processed/04_sarek_joakim \
    --genome GATK.GRCh38 \
    --igenomes_base /data/local/reference/aws/igenomes \
    --max_cpus 63 \
    --tools mutect2,controlfreec,cnvkit,msisensorpro,vep,tiddit,haplotypecaller,manta,snpeff \
    --wes


nextflow run nf-core/sarek \
    -profile singularity \
    -c custom.config \
    -w work/batch_5 \
    -resume \
    --intervals /data/local/reference/agilent/SureSelect_XT_HS_Human_All_Exon_V8/hg38/S33266340_Padded.reformatted.bed \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --max_memory '120.GB' \
    --input inputs/batch_5.csv \
    --outdir /data/local/MD_project/data/exome/processed/04_sarek_joakim \
    --genome GATK.GRCh38 \
    --igenomes_base /data/local/reference/aws/igenomes \
    --max_cpus 63 \
    --tools mutect2,controlfreec,cnvkit,msisensorpro,vep,tiddit,haplotypecaller,manta,snpeff \
    --wes
