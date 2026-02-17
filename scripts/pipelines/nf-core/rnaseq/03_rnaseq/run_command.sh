#!/bin/bash


nextflow run nf-core/rnaseq \
--input /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/helper_scripts/rnaseq_final.csv \
--outdir /home/user_oruko/work/processed/rna \
--save_merged_fastq true \
-profile singularity -r 3.18.0 \
--fasta /home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
--gtf /home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.108.gtf \
--bbsplit_fasta_list  /home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/fasta_list.csv --save_reference \
--skip_bbsplit false --save_bbsplit_reads --skip_preseq false \
-c rna.config --rseqc_modules bam_stat,inner_distance,infer_experiment,junction_saturation,read_distribution,read_duplication \
-resume 