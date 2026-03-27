#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.TUMI_reads= "/media/cph/Store4-USB/kevin/exome/01_trimmed_umis/*PDX*_R{1,2}_*.fastq.gz"


//input_disambiguate_ch = bam_human_ch.join(bam_mouse_ch, by: 0).view()

process BBSPLIT {
    containerOptions '-B /media/cph/Store4-USB/kevin/references:/media/cph/Store4-USB/kevin/references'
    maxForks 2
    publishDir params.bbsplit_dir, mode: 'copy'
    input:
        file TUMI_read_pair
    output:
        file "*"
    //publishDir "${params.outdir}/CRAM_TO_BAM_HUMAN", pattern: "*.*", mode: 'symlink'

    script:
    
    """
    prefix=\$(basename ${TUMI_read_pair[1]} | sed 's/_S.*//')
    /opt/BBMap/bbmap/bbsplit.sh  in=${TUMI_read_pair[1]} in2=${TUMI_read_pair[2]} ref_x=${params.fasta_human} ref_y=${params.fasta_mouse} out_x=\${prefix}_human_R#.fastq.gz out_y=\${prefix}_mouse.fastq.gz out_ambig=\${prefix}_ambig.fastq.gz out_unmatched=\${prefix}_unmatched.fastq.gz
    """
}


workflow {
    read_pairs_ch = Channel.fromFilePairs(params.TUMI_reads, flat: true)
    BBSPLIT(read_pairs_ch)
}