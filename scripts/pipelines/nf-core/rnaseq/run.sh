#!/usr/bin/env bash


#01_TRIM_UMIs
RAW_reads_dir="/home/user_oruko/data/raw/rna/all_fastqs/"
EXOME_SIF=/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome_2.sif
TRIM_UMIs=/home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/01_trim_umis/trim_umis.nf

TRIM_UMIS_WORK=/home/user_oruko/work/processed/rna/work/01_TRIM_UMIs
TRIMMED_DIR=/home/user_oruko/work/processed/rna/01_TRIM_UMIs
SIF=/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome_2.sif

#nextflow run $TRIM_UMIs -with-singularity $SIF --TRIMMED_DIR $TRIMMED_DIR  -work-dir $TRIM_UMIS_WORK



#RNASEQ

#nextflow run nf-core/rnaseq \
#--input /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/helper_scripts/final_rnaseq_tumi.csv \
#--outdir /home/user_oruko/work/processed/rna_trim_umis/rnaseq \
#--save_merged_fastq true \
#-profile singularity -r 3.18.0 \
#--fasta /home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
#--gtf /home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.108.gtf \
#--bbsplit_fasta_list  /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/fasta_list.csv --save_reference \
#--skip_bbsplit false --save_bbsplit_reads --skip_preseq false \
#-c /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/03_rnaseq/rna.config --rseqc_modules bam_stat,inner_distance,infer_experiment,junction_saturation,read_distribution,read_duplication \
#-resume \
#-w /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/work_rnaseq

#RNAVAR

nextflow run nf-core/rnavar -profile singularity --input /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/helper_scripts/final_rnavar_tumi.csv --outdir "/home/user_oruko/work/processed/rna_trim_umis/rnavar_nfcore" --genome GRCh38 \
--dbsnp '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf' \
--dbsnp_tbi '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx' \
--known_indels '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' \
--known_indels_tbi '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi' \
-w /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/work_rnavar