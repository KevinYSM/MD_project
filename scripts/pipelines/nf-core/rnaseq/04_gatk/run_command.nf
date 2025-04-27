#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/standard_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/"

params.RNA_raw_reads_directory='/home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/04_gatk/fastqs/*{1,2}.fastq.gz'
params.EXOME_dir='/home/user_oruko/data/raw/exome/GWA-JN-374/fastqs/*_R{1,2}*.fastq.gz'
params.fasta="/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
params.gtf="/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.108.gtf"
params.genomeDir="/home/user_oruko/other/STAR_genome_dir"


//Map to reference using STAR

def check_params() {

    if( params.remove('help') ) {
        params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
        exit 0
    }

    // additional validation here
}

process STAR_align_human{
        maxForks 3
        publishDir '/home/user_oruko/work/processed/rna/gatk/STAR_align_human_2', mode: 'copy'

        
        input:
        tuple val(sample_id), path(read_files)
        

        output:
        path '*.bam'

        script:
        """
        genomeDir=/home/user_oruko/other/STAR_genome_dir/
        sample=\$(basename "${read_files[0]}")
        ID=\$(basename ${read_files[0]})
        SM="\$sample"
        ls -l ${read_files[0]}
        ls -l ${read_files[1]}
        PL=ILLUMINA
        LB="\$SM" # Assume that each sample has only been subjected to one library preparation, which has been used for both runs (no other information has been given to indicate otherwise).
        maxReadLength=\$(zcat "${read_files[0]}" | head -40000 | awk 'NR%4==2{print length(\$0)}' | sort -nr | head -1)
        echo \$maxReadLength>>"/home/user_oruko/work/MD_project/scripts/pipelines/nf-core/rnaseq/04_gatk/sample.bam"
        
        STAR    --runThreadN 30 \
                --genomeDir "/home/user_oruko/other/STAR_genome_dir/" \
                --readFilesIn ${read_files[0]} ${read_files[1]} \
                --outFileNamePrefix \$sample \
                --outSAMtype BAM Unsorted SortedByCoordinate \
                --outSAMunmapped Within \
                --outReadsUnmapped Fastx \
                --quantMode TranscriptomeSAM GeneCounts \
                --twopassMode Basic \
                --outFilterType BySJout \
                --outSAMattrRGline "ID:\${ID} PL:\${PL} LB:\${LB} SM:\${SM}" \
                --sjdbOverhang \$(expr \$maxReadLength - 1) \
                --sjdbGTFfile "${params.gtf}" \
                --outSAMmapqUnique 60 \
                --readFilesCommand zcat
        """

}

//Generate unaligned bam files
process UnalignedBams{
        maxForks 5
        publishDir '/home/user_oruko/work/processed/rna/gatk/BAM_unaligned', mode: 'copy'

        
        input:
        tuple val(sample_id), path(read_files)
        

        output:
        path '*.bam'
        script:
        """
        #Patient IDs with read group name S152
        declare -a patient_ids_152=("9" "16" "20" "28" "30" "31" "36" "39" "50" "54")

        
        patient_id=\$(echo "${read_files[0]}" | sed 's/p\\([0-9]\\{1,3\\}\\)_primary.*\$/\\1/')
        sample_group="S301"    
        for id in \${patient_ids_152[@]}; do
            
            if [[ \${patient_id} == "\${id}" ]]; then
                sample_group="S152"
                
            fi
        done


        java -jar /home/user_oruko/other/picard/picard.jar FastqToSam \
        F1=${read_files[0]} \
        F2=${read_files[1]} \
        O=unaligned_"${sample_id}".bam \
        SM=\${patient_id} \
        READ_GROUP_NAME=\${sample_group} \
        PLATFORM=ILLUMINA
        
        """
}
//Data cleanup using MergeBamAlignment and Markduplicates
process MergeBamAlignment{
    input:
        file BAM_file
        

    output:
        path "*.bam"

        """
        base=\$(basename $BAM_file)
        java -jar /picard/picard.jar MergeBamAlignment \
        --ALIGNED_BAM $BAM_file \\
        UNMAPPED=unmapped.bam \\
        O=merge_alignments.bam \\
        R=reference_sequence.fasta

        
        """
    
}
//SplitNCigarReads

//Base Quality Recalibration

//Variant Calling

//Variant Filtering


workflow {

    RNA_RAW_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, suffix: '.fastq.gz')
    RNA_RAW_ch.view()  // To check the files are properly captured by the channel
    //STAR_align_human(RNA_RAW_ch)
    UnalignedBams(RNA_RAW_ch)
}
