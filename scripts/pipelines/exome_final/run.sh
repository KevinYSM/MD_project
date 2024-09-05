#!/bin/bash
BASE_DIR="/data/local/MD_project"
JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/21.0.3-tem"
JAVA_CMD="/home/ubuntu/.sdkman/candidates/java/21.0.3-tem"


SIF="/data/local/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome_2.sif"
EXOME_RAW_READS="/data/local/MD_project/data/exome/raw/*/fastqs/*_R{1,2}_*.fastq.gz"
TRIMMED_DIR="/data/local/MD_project/data/exome/processed_final/01_trimmed_umis_redone_"



#01_Trim UMIs
TRIM_UMIS_NF="/data/local/MD_project/scripts/pipelines/exome_final/steps/01_trim_umis/trim_umis.nf"
TRIM_UMIS_WORK="/home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/work/01_trim_umis"
#nextflow run $TRIM_UMIS_NF -with-singularity $SIF -resume  --TRIMMED_DIR $TRIMMED_DIR  -work-dir $TRIM_UMIS_WORK


#Generate Tumour Sample Sheets
#python3 /data/local/MD_project/scripts/pipelines/exome_final/helper_scripts/get_tumi_tumour.py $TRIMMED_DIR

#02_Disambiguate
WORK_HUMAN_DIR="/data/local/MD_project/scripts/pipelines/exome_final/work/02_disambiguate/human"
RESULTS_HUMAN_DIR="/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_human"

WORK_MOUSE_DIR="/data/local/MD_project/scripts/pipelines/exome_final/work/02_disambiguate/mouse"
RESULTS_MOUSE_DIR="/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_mouse"
RESULTS_DISAMBIGUATE_DIR="/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_disambiguate"

MAP_HUMAN="/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/map_human.sh"
MAP_MOUSE="/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/map_mouse.sh"

#Please note, this step below only runs on tumour samples:
#source $MAP_HUMAN $WORK_HUMAN_DIR $RESULTS_HUMAN_DIR
#source $MAP_MOUSE $WORK_MOUSE_DIR $RESULTS_MOUSE_DIR

CRAM_HUMAN='/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_human/preprocessing/mapped/*/*.cram'
CRAM_MOUSE='/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_mouse/preprocessing/mapped/*/*.cram'
FASTA_MOUSE="/data/local/reference/aws/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
FASTA_HUMAN="/data/local/reference/aws/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"


DISAMBIGUATE_NF="/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/disambiguate.nf"



#nextflow run $DISAMBIGUATE_NF -with-singularity $SIF --cram_human "$CRAM_HUMAN" --cram_mouse "$CRAM_MOUSE" \
# --fasta_mouse $FASTA_MOUSE \
# --fasta_human $FASTA_HUMAN \
# --outdir $RESULTS_DISAMBIGUATE_DIR -resume \
# --max_memory '120.GB' \
# --max_cpus 63 \


#03_Prepare Samplesheets
samplesheet_path="/data/local/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/exome_final_samplesheet.csv"
#SAREK
#nextflow run nf-core/sarek -with-singularity $SIF --input $samplesheet_path \
#-resume \
# -profile singularity --fasta "/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
#--tools mutect2,haplotypecaller,ascat,cnvkit,msisensorpro --wes \
#--outdir "/data/local/MD_project/data/exome/processed_final/sarek" 

batch_1="/data/local/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/batches/batch_1.csv"
batch_2="/data/local/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/batches/batch_2.csv"
batch_4="/data/local/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/batches/batch_4.csv"
#nextflow run nf-core/sarek -with-singularity $SIF --input $batch_4 \
#-resume \
# -profile singularity --fasta "/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
# -work-dir "/data/local/MD_project/scripts/pipelines/exome_final/work/batch_4" \
#--tools mutect2,haplotypecaller,ascat,cnvkit,msisensorpro --wes \
#--outdir "/data/local/MD_project/data/exome/processed_final/sarek" 

batch_base="/data/local/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/batches"
work_dir_base="/data/local/MD_project/scripts/pipelines/exome_final/work"
out_dir="/data/local/MD_project/data/exome/processed_final/sarek"
fasta_file="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

# Loop through batch_1 to batch_8
#for i in {2..8}; do
#    batch_file="${batch_base}/batch_${i}.csv"
#    work_dir="${work_dir_base}/batch_${i}"
#  
#    nextflow run nf-core/sarek -with-singularity $SIF --input $batch_file \
#    -profile singularity --fasta $fasta_file \
#    -work-dir $work_dir \
#    --tools mutect2,haplotypecaller,ascat,cnvkit,msisensorpro --wes \
#    --outdir $out_dir
#done

#VCF
MUTECT2="/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2/*/*filtered.vcf.gz"
export NXF_SINGULARITY_OPTS="--bind /data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2:/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2"
export SINGULARITY_BINDPATH="/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2:/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2"
export SINGULARITY_BIND="/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2:/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2"
if [ ! -d output_vep_updated ]
    then
        mkdir output_vep_updated
fi

#01_gatk_filter_mutect_calls
#nextflow run /data/local/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/01_gatk_filter_mutect_calls.nf -with-singularity $SIF \
# --max_memory '120.GB' \
# --max_cpus 63 -work-dir "/data/local/MD_project/scripts/pipelines/exome_final/work/gatk_filter"

 #02_vep
nextflow run /data/local/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/02_vep.nf  -with-singularity $SIF  \
--max_memory '120.GB' \
--max_cpus 63 -work-dir "/data/local/MD_project/scripts/pipelines/exome_final/work/vep" 

#03_vcf2maf
nextflow run /data/local/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/03_vcf2maf.nf  -with-singularity $SIF  \
--max_memory '120.GB' \
 --max_cpus 63 \
-c /data/local/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf.config -work-dir "/data/local/MD_project/scripts/pipelines/exome_final/work/vcf2maf"

