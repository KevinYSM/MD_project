#!/usr/bin/env nextflow
params.unfiltered_vcfs="/home/user_oruko/work/processed/exome/vcf_processing_intervals/output_vep_refseq/*.ann.vcf"
params.outdir="/home/user_oruko/work/processed/exome/vcf_processing_intervals/filtered_vcfs_new"


process filterPASS{
    input:
        path unfiltered_vcfs
    output:
        path "*.PASS.filtered.vcf", emit: PASS_vcf
    publishDir "/home/user_oruko/work/processed/exome/vcf_processing_intervals/filtered_vcfs_PASS_2", mode: 'copy'
    script:

    """
    basename=\$(basename "${unfiltered_vcfs}" ".vcf")
    vcftools --vcf ${unfiltered_vcfs} --remove-filtered-all --recode-INFO-all --recode --stdout > "\$basename".PASS.filtered.vcf
    """
}
process filterDepth {
    containerOptions "--bind /home/user_oruko/other/:/home/user_oruko/other/"
    input:
        path PASS_vcfs
    
    script:
    """
    basename=\$(basename "${PASS_vcfs}" ".vcf")
    java -jar /home/user_oruko/other/vcffilter-assembly-0.2.jar  -I ${PASS_vcfs} -o /home/user_oruko/work/processed/exome/vcf_processing_intervals/filtered_vcfs_new/"\$basename".filtered.vcf --minSampleDepth 10 --minTotalDepth 30 
    """
}


workflow {
    VCF_ch = Channel.fromPath(params.unfiltered_vcfs, checkIfExists: true)
    filterPASS(VCF_ch)
    //filterDepth(filterPASS(VCF_ch))
}
