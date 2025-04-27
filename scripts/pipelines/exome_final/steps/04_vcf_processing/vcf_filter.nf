params.unfiltered_vcfs="/home/user_oruko/work/processed/exome/vcf_processing_intervals/output_vep_refseq/*.ann.vcf"
params.outdir="/home/user_oruko/work/processed/exome/vcf_processing_intervals/filtered_vcfs_qual_30"

process VCF2MAF {
    containerOptions "--bind /home/user_oruko/other/:/home/user_oruko/other/"
    input:
        path unfiltered_vcfs
    publishDir "${params.outdir}", mode: 'copy'
    script:
    """
    basename=\$(basename "${unfiltered_vcfs}" ".vcf")
    java -jar /home/user_oruko/other/vcffilter-assembly-0.2.jar  -I ${unfiltered_vcfs} -o /home/user_oruko/work/processed/exome/vcf_processing_intervals/filtered_vcfs_qual_30/"\$basename".filtered.vcf --minSampleDepth 10 --minTotalDepth 30 
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
