#!/bin/bash
BASE_DIR="/home/kevin/MD_project"
DATA_DIR="/media/cph/Store4-USB/kevin/exome"
JAVA_HOME="/home/kevin/bin/.sdkman/candidates/java/25.0.2-tem"
JAVA_CMD="/home/kevin/bin/.sdkman/candidates/java/25.0.2-tem"


#singularity build "${DATA_DIR}/00_prep/exome.sif" "${BASE_DIR}/scripts/pipelines/exome_final/steps/00_prep/singularity/exome.def"
#singularity build "${DATA_DIR}/00_prep/exome_vep.sif" "${BASE_DIR}/scripts/pipelines/exome_final/steps/00_prep/singularity/exome_vep.def"

SIF="${DATA_DIR}/00_prep/exome_vep.sif"
EXOME_RAW_READS="${DATA_DIR}/data/raw/exome/*/fastqs/*_R{1,2}_*.fastq.gz"
TRIMMED_DIR="${DATA_DIR}/01_trimmed_umis"



#01_Trim UMIs
TRIM_UMIS_NF="${BASE_DIR}/scripts/pipelines/exome_final/steps/01_trim_umis/trim_umis.nf"
TRIM_UMIS_WORK="${DATA_DIR}/work/01_trim_umis"
#nextflow run $TRIM_UMIS_NF -with-singularity $SIF -resume  --TRIMMED_DIR $TRIMMED_DIR  -work-dir $TRIM_UMIS_WORK 

#Generate Tumour Sample Sheets
#python3 ${BASE_DIR}/scripts/pipelines/exome_final/helper_scripts/get_tumi_tumour.py $TRIMMED_DIR

#02_Disambiguate
# Two parts: mapping to human and mouse genomes, and then running disambiguate on the resulting CRAM files. The first part is run on tumour samples only, while the second part is run on all samples (tumour and normal).
WORK_HUMAN_DIR="${DATA_DIR}/work/02_disambiguate/human"
RESULTS_HUMAN_DIR="${DATA_DIR}/02_disambiguate/results_human"

WORK_MOUSE_DIR="${DATA_DIR}/work/02_disambiguate/mouse"
RESULTS_MOUSE_DIR="${DATA_DIR}/02_disambiguate/results_mouse"
RESULTS_BBSPLIT_DIR="${DATA_DIR}/02_bbsplit/results_bbsplit"

MAP_HUMAN="${BASE_DIR}/scripts/pipelines/exome_final/steps/02_disambiguate/map_human.sh"
MAP_MOUSE="${BASE_DIR}/scripts/pipelines/exome_final/steps/02_disambiguate/map_mouse.sh"



FASTA_MOUSE="/media/cph/Store4-USB/kevin/references/mouse/ensembl/Mus_musculus.GRCm39.dna.primary_assembly.fa"
FASTA_HUMAN="/media/cph/Store4-USB/kevin/references/human/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

BBSPLIT_NF_1="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02_bbsplit/bbsplit_50.nf"
BBSPLIT_NF_2="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02_bbsplit/bbsplit_73.nf"

#nextflow run $BBSPLIT_NF_1 -with-singularity $SIF \
#--fasta_human $FASTA_HUMAN \
#--fasta_mouse $FASTA_MOUSE \
#--bbsplit_dir $RESULTS_BBSPLIT_DIR -resume \
#--max_memory '100.GB' \
#-w /media/cph/Store4-USB/kevin/exome/02_bbsplit/work



#03_SAREK
SAREK_SAMPLESHEET_1="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/sarek_1.csv"
SAREK_SAMPLESHEET_2="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/sarek_2.csv"
SAREK_SAMPLESHEET_3="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/sarek_3.csv"
SAREK_SAMPLESHEET_4="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/sarek_4.csv"
SAREK_SAMPLESHEET_5="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/sarek_5.csv"
SAREK_SAMPLESHEET_50_73="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/sarek_50_73.csv"


#nextflow run nf-core/sarek -with-singularity $SIF --input $SAREK_SAMPLESHEET_50_73 \
#    -w "/media/cph/Store4-USB/kevin/exome/03_sarek/work" \
#    -profile singularity \
#    --intervals "/media/cph/Store4-USB/kevin/references/human/agilent/S33266340_hg38/S33266340_Padded.sorted.merged.bed " \
#    --cf_contamination_adjustment FALSE \
#    --cf_contamination 0 \
#    --genome GATK.GRCh38 \
#    --tools mutect2,vep,haplotypecaller,cnvkit \
#    --wes \
#    --outdir "/media/cph/Store4-USB/kevin/exome/03_sarek" \
#    -c /home/kevin/MD_project/scripts/pipelines/exome_final/steps/03_sarek/sarek.config \
#    -resume 




#vcf processing
MUTECT2="/media/cph/Store4-USB/kevin/exome/03_sarek/variant_calling/mutect2/*/*filtered.vcf.gz"
export NXF_SINGULARITY_OPTS="--bind /media/cph/Store4-USB/kevin/exome/03_sarek/variant_calling/mutect2/:/media/cph/Store4-USB/kevin/exome/03_sarek/variant_calling/mutect2/"
export SINGULARITY_BINDPATH="/media/cph/Store4-USB/kevin/exome/03_sarek/variant_calling/mutect2:/media/cph/Store4-USB/kevin/exome/03_sarek/variant_calling/mutect2"
export SINGULARITY_BIND="/media/cph/Store4-USB/kevin/exome/03_sarek/variant_calling/mutect2:/media/cph/Store4-USB/kevin/exome/03_sarek/variant_calling/mutect2"
if [ ! -d output_vep_updated ]
    then
        mkdir output_vep_updated
fi



nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vep.nf  -with-singularity $SIF  \
--max_memory '185.GB' \
 --max_cpus 94 \
 -w "/media/cph/Store4-USB/kevin/exome/04_vcf_processing/vep/work " 

nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf_filter.nf \
 --max_memory '185.GB' --max_cpus 94 \



nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/add_af.nf \
--max_memory '185.GB' \
--max_cpus 94 

nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf2maf.nf  -with-singularity $SIF  \
--max_memory '185.GB' \
--max_cpus 94 \
-c /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf.config