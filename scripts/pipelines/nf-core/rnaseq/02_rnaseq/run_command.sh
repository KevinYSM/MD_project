#!/bin/bash
nextflow run nf-core/rnaseq \
--input /data/local/MD_scholarly/scripts/pipelines/nf-core/rnaseq/helper_scripts/19_samples.csv \
--outdir /data/local/MD_scholarly/data/processed/nf-core/rnaseq 
--save_merged_fastq true \
-profile singularity -r 3.12.0 \ 
--bbsplit_fasta_list fasta_list.csv \
 —-bbsplit_fasta_list  /data/local/MD_scholarly/scripts/pipelines/nf-core/rnaseq/fasta_list.csv --save_reference\
 --skip_bbsplit false --save_bbsplit_reads --skip_preseq false \
--fasta /data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
--gtf /data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf \
-c rna.config --rseqc_modules bam_stat,inner_distance,infer_experiment,junction_saturation,read_distribution,read_duplication \
-resume \
/