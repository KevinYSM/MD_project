



#02b_Prepare Samplesheets
samplesheet_path="/home/kevin/MD_project/scripts/pipelines/exome_final/steps/02b_samplesheets/samplesheets/exome_final_MASTER.csv"

#03_SAREK
nextflow run nf-core/sarek -with-singularity $SIF --input $samplesheet_path \
    -w "/home/kevin/MD_project/scripts/pipelines/exome_final/work_intervals/sarek" \
    -profile singularity \
    --intervals "/home/kevin/MD_project/data/references/agilent/S33266340_Padded_merge.clean.bed" \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --genome GATK.GRCh38 \
    --tools mutect2,msisensorpro,vep,haplotypecaller, cnvkit \
    --wes \
    --outdir "/home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek_intervals" \
    - resume \
    -c "/home/kevin/MD_project/scripts/pipelines/exome_final/configs/sarek.config" 



#vcf processing
MUTECT2="/home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek/variant_calling/mutect2/*/*filtered.vcf.gz"
export NXF_SINGULARITY_OPTS="--bind /home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek/variant_calling/mutect2:/home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek/variant_calling/mutect2"
export SINGULARITY_BINDPATH="/home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek/variant_calling/mutect2:/home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek/variant_calling/mutect2"
export SINGULARITY_BIND="/home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek/variant_calling/mutect2:/home/kevin/MD_project/scripts/pipelines/exome_final/processed/exome/sarek/variant_calling/mutect2"
if [ ! -d output_vep_updated ]
    then
        mkdir output_vep_updated
fi



nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vep.nf  -with-singularity $SIF  \
--max_memory '185.GB' \
 --max_cpus 94 \
 -w "/home/kevin/MD_project/scripts/pipelines/exome_final/work_intervals/vep_refseq" 

nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf_filter.nf \
 --max_memory '185.GB' --max_cpus 94 \



nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/add_af.nf \
--max_memory '185.GB' \
--max_cpus 94 

nextflow run /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf2maf.nf  -with-singularity $SIF  \
--max_memory '185.GB' \
--max_cpus 94 \
-c /home/kevin/MD_project/scripts/pipelines/exome_final/steps/04_vcf_processing/vcf.config