#!/bin/bash
BASE_DIR="/home/user_oruko/work/MD_project"
JAVA_HOME="/home/user_oruko/other/bin/.sdkman/candidates/java/21.0.5-tem"
JAVA_CMD="/home/user_oruko/other/bin/.sdkman/candidates/java/21.0.5-tem"


SIF="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome_2.sif"
EXOME_RAW_READS="/home/user_oruko/data/raw/exome/*/fastqs/*_R{1,2}_*.fastq.gz"
TRIMMED_DIR="/home/user_oruko/work/processed/exome/01_trimmed_umis"



#01_Trim UMIs
TRIM_UMIS_NF="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/01_trim_umis/trim_umis.nf"
TRIM_UMIS_WORK="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/work/01_trim_umis_re_demux"
nextflow run $TRIM_UMIS_NF -with-singularity $SIF -resume  --TRIMMED_DIR $TRIMMED_DIR  -work-dir $TRIM_UMIS_WORK


#Generate Tumour Sample Sheets
python3 /home/user_oruko/work/MD_project/scripts/pipelines/exome_final/helper_scripts/get_tumi_tumour.py $TRIMMED_DIR

#02_Disambiguate
# Two parts: mapping to human and mouse genomes, and then running disambiguate on the resulting CRAM files. The first part is run on tumour samples only, while the second part is run on all samples (tumour and normal).
WORK_HUMAN_DIR="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/work/02_disambiguate/human"
RESULTS_HUMAN_DIR="/home/user_oruko/work/processed/exome/02_disambiguate/results_human"

WORK_MOUSE_DIR="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/work/02_disambiguate/mouse"
RESULTS_MOUSE_DIR="/home/user_oruko/work/processed/exome/02_disambiguate/results_mouse"
RESULTS_DISAMBIGUATE_DIR="/home/user_oruko/work/processed/exome/02_disambiguate/results_disambiguate"

MAP_HUMAN="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/map_human.sh"
MAP_MOUSE="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/map_mouse.sh"

#Please note, this step below only runs on tumour samples:
source $MAP_HUMAN $WORK_HUMAN_DIR $RESULTS_HUMAN_DIR
source $MAP_MOUSE $WORK_MOUSE_DIR $RESULTS_MOUSE_DIR

CRAM_HUMAN='/home/user_oruko/work/processed/exome/02_disambiguate/results_human/preprocessing/mapped/*/*.cram'
CRAM_MOUSE='/home/user_oruko/work/processed/exome/02_disambiguate/results_mouse/preprocessing/mapped/*/*.cram'


FASTA_MOUSE="/home/user_oruko/data/references/aws/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
FASTA_HUMAN="/home/user_oruko/data/references/aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"


DISAMBIGUATE_NF="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/disambiguate.nf"

nextflow run $DISAMBIGUATE_NF -with-singularity $SIF --cram_human "$CRAM_HUMAN" --cram_mouse "$CRAM_MOUSE" \
 --fasta_mouse $FASTA_MOUSE \
 --fasta_human $FASTA_HUMAN \
 --outdir $RESULTS_DISAMBIGUATE_DIR -resume \
 --max_memory '175.GB' \
 --max_cpus 93 \


#02b_Prepare Samplesheets
samplesheet_path="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/exome_final_MASTER.csv"

#03_SAREK
nextflow run nf-core/sarek -with-singularity $SIF --input $samplesheet_path \
    -w "/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/work_intervals/sarek" \
    -profile singularity \
    --intervals "/home/user_oruko/data/references/agilent/S33266340_Padded_merge.clean.bed" \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --genome GATK.GRCh38 \
    --tools mutect2,msisensorpro,vep,haplotypecaller, cnvkit \
    --wes \
    --outdir "/home/user_oruko/work/processed/exome/sarek_intervals" \
    - resume \
    -c "/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/configs/sarek.config" 



#vcf processing
MUTECT2="/home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2/*/*filtered.vcf.gz"
export NXF_SINGULARITY_OPTS="--bind /home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2:/home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2"
export SINGULARITY_BINDPATH="/home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2:/home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2"
export SINGULARITY_BIND="/home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2:/home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2"
if [ ! -d output_vep_updated ]
    then
        mkdir output_vep_updated
fi



nextflow run /home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vep.nf  -with-singularity $SIF  \
--max_memory '185.GB' \
 --max_cpus 94 \
 -w "/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/work_intervals/vep_refseq" 

nextflow run /home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf_filter.nf \
 --max_memory '185.GB' --max_cpus 94 \



nextflow run /home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/add_af.nf \
--max_memory '185.GB' \
--max_cpus 94 

nextflow run /home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf2maf.nf  -with-singularity $SIF  \
--max_memory '185.GB' \
--max_cpus 94 \
-c /home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf.config