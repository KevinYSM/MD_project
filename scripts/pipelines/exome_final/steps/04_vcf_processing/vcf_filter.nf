params.unfiltered_vcfs="/data/local/MD_project/data/output_vep/*.ann.vcf"
params.outdir="/data/local/MD_project/data/filtered_vcfs"

process VCF2MAF {
    input:
        path unfiltered_vcfs
    publishDir "${params.outdir}", mode: 'copy'
    script:
    """
    basename=\$(basename "${unfiltered_vcfs}" ".vcf")
    vcffilter -f "(DP > 5 & MMQ > 30 & )" ${unfiltered_vcfs} > /data/local/MD_project/data/filtered_vcfs/"\$basename".filtered.vcf
    """
}


workflow {
    VCF_ch = Channel.fromPath(params.unfiltered_vcfs, checkIfExists: true)
    
    view(VCF_ch)  // To check the files are properly captured by the channel
    
    VCF_ch.subscribe { file -> 
        println "Processing file: $file"
    }

    VCF2MAF(VCF_ch)
}
