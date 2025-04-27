#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/standard_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/"

params.RNA_raw_reads_directory='/home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/04_gatk/fastqs/*{1,2}.fastq.gz'
params.EXOME_dir='/home/user_oruko/data/raw/exome/GWA-JN-374/fastqs/*_R{1,2}*.fastq.gz'
params.fasta="/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
params.gtf="/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.108.gtf"
params.genomeDir="/home/user_oruko/other/STAR_genome_dir"
params.RNA_unaligned_bams_directory="/home/user_oruko/work/processed/rna/gatk/BAM_unaligned/*primary.bam"
params.RNA_aligned_bams_directory="/home/user_oruko/work/processed/rna/gatk/BAM_aligned/*Aligned.out.bam"

params.RNA_merge_bam_aligned="/home/user_oruko/work/processed/rna/gatk/BAM_replaced_read_groups_2/broken_rg/output.bam"


//Data cleanup using MergeBamAlignment and Markduplicates
process MergeBamAlignment{
    publishDir "/home/user_oruko/work/processed/rna/gatk/BAM_merge_alignment", mode: 'copy'
    input:
        file aligned_BAM_file

    output:
        path '*.bam'
    script:
        """
        
        filename=\$(basename "$aligned_BAM_file")
        patient_id=\${filename%%_p*}
        
        start_term="unaligned_"
        end_term="_primary"
       
        for f in ${params.RNA_unaligned_bams_directory}; do

            if [[ "\$f" =~ "\${start_term}\${patient_id}\${end_term}" ]]; then
                

                java -Djava.io.tmpdir=/home/user_oruko/work/tmp -jar /home/user_oruko/other/picard/picard.jar MergeBamAlignment \
                --ALIGNED_BAM $aligned_BAM_file \\
                --UNMAPPED_BAM "\$f" \\
                --OUTPUT "\$patient_id"_merge_alignment.bam \\
                --REFERENCE_SEQUENCE ${params.fasta}
            fi
        done
        """
}

process MarkDuplicates {
    maxForks 5
    publishDir "/home/user_oruko/work/processed/rna/gatk/BAM_duplicates", mode: 'copy'
    input:
        file ALIGNED_BAM
    output:
        path ('*')
    script:
    """

    filename=\$(basename "${ALIGNED_BAM}")

    java -Djava.io.tmpdir=/home/user_oruko/work/tmp -jar /home/user_oruko/other/picard/picard.jar MarkDuplicates \
      I=${ALIGNED_BAM} \
      O=\$filename.duplicates.bam \
      M=\$filename.marked_dup_metrics.txt
    """
}

process AddReadGroups {
    publishDir "/home/user_oruko/work/processed/rna/gatk/BAM_replaced_read_groups", mode: 'copy'
    input:
        file ALIGNED_BAM
    output:
        path ('*')
    script:
    """

    filename=\$(basename "${ALIGNED_BAM}")
    #Patient IDs with read group name S152
    declare -a patient_ids_152=("9" "16" "20" "28" "30" "31" "36" "39" "50" "54")
    
    
    patient_id=\$(echo "${ALIGNED_BAM}" | sed 's/p\\([0-9]\\{1,3\\}\\)_primary.*\$/\\1/')
    sample_group="S301"    
    for id in \${patient_ids_152[@]}; do
        
        if [[ \${patient_id} == "\${id}" ]]; then
            sample_group="S152"
            
        fi
    done
    SM="\${patient_id}""_\${sample_group}"
    java -Djava.io.tmpdir=/home/user_oruko/work/tmp -jar /home/user_oruko/other/picard/picard.jar AddOrReplaceReadGroups \
      I=${ALIGNED_BAM} \
      O=\$filename.rg.bam \
      RGID=\$sample_group \
      RGLB="unchecked" \
      RGPL='illumina' \
      RGPU="unchecked" \
      RGSM=\$SM
    """

}
workflow {

    //unaligned_bams_ch=Channel.fromPath(params.RNA_unaligned_bams_directory)

    //aligned_bams_ch=Channel.fromPath(params.RNA_aligned_bams_directory)

    //MergeBamAlignment(aligned_bams_ch)
    RNA_merged_bam_ch=Channel.fromPath(params.RNA_merge_bam_aligned)
    //AddReadGroups(RNA_merged_bam_ch)
    MarkDuplicates(RNA_merged_bam_ch)
    //unaligned_bams_ch.view()
    //aligned_bams_ch.view()

}
