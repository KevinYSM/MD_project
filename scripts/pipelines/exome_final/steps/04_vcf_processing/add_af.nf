#!/usr/bin/env nextflow
params.pass_vcfs="/home/user_oruko/work/processed/exome/vcf_processing_intervals/filtered_vcfs_PASS_2/*.vcf"
process addAF{
    input:
        path PASS_vcfs
    output:
        path "*"
    publishDir "/home/user_oruko/work/processed/exome/vcf_processing_intervals/filtered_vcf_PASS_AF", mode: 'copy'
    script:

    """
    basename=\$(basename "${PASS_vcfs}" ".vcf")
    bcftools +fill-tags ${PASS_vcfs} -o "\$basename".AF.vcf -- -t AF  
    """
}


workflow {
    VCF_ch = Channel.fromPath(params.pass_vcfs, checkIfExists: true)
    addAF(VCF_ch)
    //filterDepth(filterPASS(VCF_ch))
}
