#!/usr/bin/env nextflow

nextflow.enable.dsl=2


cram_human_ch = Channel.fromFilePairs( params.cram_human, size: 1)
cram_mouse_ch = Channel.fromFilePairs( params.cram_mouse, size: 1)

//input_disambiguate_ch = bam_human_ch.join(bam_mouse_ch, by: 0).view()

process CRAM_TO_BAM_HUMAN {

    //publishDir "${params.outdir}/CRAM_TO_BAM_HUMAN", pattern: "*.*", mode: 'symlink'

    input:
        tuple val(sname), path(cram)

    output:
        tuple val(sname), path("*.bam"), emit: cram_to_bam_human_ch

    script:
    """
    samtools view -b -o ${sname}.human.bam -T ${params.fasta_human} ${cram}
    """
}

process CRAM_TO_BAM_MOUSE {

    //publishDir "${params.outdir}/CRAM_TO_BAM_MOUSE", pattern: "*.*", mode: 'symlink'

    input:
        tuple val(sname), path(cram)

    output:
        tuple val(sname), path("*.bam"), emit: cram_to_bam_mouse_ch

    script:
    """
    samtools view -b -o ${sname}.mouse.bam -T ${params.fasta_mouse} ${cram}
    """
}

process DISAMBIGUATE {

    publishDir "${params.outdir}/disambiguate", pattern: "*.*", mode: 'symlink'

    input:
        tuple val(sname), path(bam_human), path(bam_mouse)

    output:
        tuple val(sname), path("*.disambiguatedSpeciesA.bam"), emit: bam_human_disambiguated_ch

    script:
    """
    ngs_disambiguate -s "${sname}" -o "./" -a bwa "${bam_human}" "${bam_mouse}"
    """
}

process BAM_TO_FASTQ {
    
    publishDir "${params.outdir}/bam_to_fastq", pattern: "*.*", mode: 'move'

    input:
        tuple val(sname), path(bam)

    output:
        tuple val(sname), path("${sname}_human_R1.fastq.gz"), path("${sname}_human_R2.fastq.gz"), emit:bam_to_fastq_ch

    script:
    """
    TMP_DIR=./\$(basename "${bam}" ".bam")_tmp

    # Sort by read name
    samtools sort \
        ${bam} \
        -o queryname_sorted.bam \
        -n
     
    #gatk SortSam \
    #    --TMP_DIR=\$TMP_DIR \
    #    -I ${bam} \
    #    -O queryname_sorted.bam \
    #    -SO queryname

    TMP_DIR_2=./\$(basename "${bam}" ".bam")_tmp_2

    # Convert to FASTQ while keeping the full reads (including soft trimmed parts)
    gatk SamToFastq \
        --TMP_DIR=\$TMP_DIR_2 \
        -I queryname_sorted.bam \
        -F ${sname}_human_R1.fastq.gz \
        -F2 ${sname}_human_R2.fastq.gz \
        --VALIDATION_STRINGENCY LENIENT \
        --INCLUDE_NON_PF_READS true

    rm queryname_sorted.bam
    """
}

workflow {


    cram_to_bam_human_ch = CRAM_TO_BAM_HUMAN(
        cram_human_ch
    )
    
    cram_to_bam_mouse_ch = CRAM_TO_BAM_MOUSE(
        cram_mouse_ch
    )

    cram_to_bam_ch = cram_to_bam_human_ch.join(cram_to_bam_mouse_ch, by: 0)
    //cram_to_bam_ch.view()

    disambiguate_ch = DISAMBIGUATE(
        cram_to_bam_ch
    )

    BAM_TO_FASTQ(
        disambiguate_ch
    )

}