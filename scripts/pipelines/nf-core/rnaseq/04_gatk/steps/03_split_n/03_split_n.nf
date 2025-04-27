mark_duplicates_loc="/home/user_oruko/work/processed/rna/gatk/BAM_duplicates/fixed_lb/*.bam"
params.fasta="/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
process SPLIT_N_TRIM {
    publishDir "/home/user_oruko/work/processed/rna/gatk/SPLIT_N", mode: 'copy'
    maxForks 6

    input:
        path(bam)
    output:
        path("*.split.bam"), emit: split_ch

    script:
    """
    gatk SplitNCigarReads  -R "${params.fasta}" -I "${bam}" -O \$(basename "${bam}" ".bam").split.bam

    """
}


workflow{
    mark_duplicates_ch=Channel.fromPath(mark_duplicates_loc)
    SPLIT_N_TRIM(mark_duplicates_ch)

}
