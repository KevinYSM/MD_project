#!/bin/bash
BASE_DIR="/data/local/MD_project"
JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/21.0.3-tem"
JAVA_CMD="/home/ubuntu/.sdkman/candidates/java/21.0.3-tem"


SIF="/data/local/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome.sif"
EXOME_RAW_READS="/data/local/MD_project/data/exome/raw/*/fastqs/*_R{1,2}_*.fastq.gz"
TRIMMED_DIR="/data/local/MD_project/data/exome/processed_final/01_trimmed_umis_redone"



#01_Trim UMIs
TRIM_UMIS_NF="/data/local/MD_project/scripts/pipelines/exome_final/steps/01_trim_umis/trim_umis.nf"
TRIM_UMIS_WORK="/home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/work/01_trim_umis"
#nextflow run $TRIM_UMIS_NF -with-singularity $SIF -resume  --TRIMMED_DIR $TRIMMED_DIR  -work-dir $TRIM_UMIS_WORK


#Generate Sample Sheets
python3 /data/local/MD_project/scripts/pipelines/exome_final/helper_scripts/04_generate_samplesheets_disambiguate.py $TRIMMED_DIR

#02_Disambiguate
WORK_HUMAN_DIR="/data/local/MD_project/scripts/pipelines/exome_final/work/02_disambiguate/human"
RESULTS_HUMAN_DIR="/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_human"

WORK_MOUSE_DIR="/data/local/MD_project/scripts/pipelines/exome_final/work/02_disambiguate/mouse"
RESULTS_MOUSE_DIR="/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_mouse"
RESULTS_DISAMBIGUATE_DIR="/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_disambiguate"

MAP_HUMAN="/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/map_human.sh"
MAP_MOUSE="/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/map_mouse.sh"
source $MAP_HUMAN $WORK_HUMAN_DIR $RESULTS_HUMAN_DIR
source $MAP_MOUSE $WORK_MOUSE_DIR $RESULTS_MOUSE_DIR

CRAM_HUMAN=$RESULTS_HUMAN_DIR+'/preprocessing/mapped/*/*.cram'
CRAM_MOUSE=$RESULTS_MOUSE_DIR+'/preprocessing/mapped/*/*.cram'

FASTA_HUMAN = '/data/local/reference/aws/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'
FASTA_MOUSE= '/data/local/reference/aws/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa'


DISAMBIGUATE_NF="/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/disambiguate.nf"



#nextflow run $DISAMBIGUATE_NF -with-singularity $SIF --cram_human $CRAM_HUMAN --cram_mouse $CRAM_MOUSE --fasta_human $FASTA_HUMAN --fasta_mouse $FASTA_MOUSE --outdir $RESULTS_DISAMBIGUATE_DIR


#03_Prepare Samplesheets

#SAREK

#vcf processing