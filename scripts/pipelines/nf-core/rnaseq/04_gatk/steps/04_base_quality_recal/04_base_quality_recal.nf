split_n_bams="/home/user_oruko/work/processed/rna/gatk/SPLIT_N/*.bam"
params.fasta="/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

params.dbsnp = '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf'
params.known_indels = '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
process baseRecalibrator{
    publishDir "/home/user_oruko/work/processed/rna/gatk/SPLIT_N", mode: 'copy'
    maxForks 6

    input:
        path(bam)
    output:
        path("*.split.bam"), emit: split_ch

    script:
    """
    gatk BaseRecalibrator \
    -I ${bam} \
    -R ${params.fasta} \
    --known-sites ${params.dbsnp} \
    --known-sites ${params.known_indels} \
    -O "\$(basename ${bam})".recal_pass1.table

    gatk ApplyBQSR \
        -I "${bam}" \
        -R "${params.fasta}" \
        --bqsr-recal-file "\$(basename ${bam})".recal_pass1.table \
        -O \$(basename "${bam}").recal.pass1.bam
    """
}

process applyRecalibration{
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

process analyzeCovariates{
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
    bams=Channel.fromPath(split_n_bams)
    base_recal_ch=baseRecalibrator(bams)
    analyze_covariates_ch=applyRecalibration(base_recal_ch)
    analyzeCovariates(analyze_covariates_ch)
}