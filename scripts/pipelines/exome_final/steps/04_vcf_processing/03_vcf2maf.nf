params.vep_files="/data/local/MD_project/data/exome/processed_final/vcf_processing_all/02_VEP/*.ann.vcf"
params.outdir="/data/local/MD_project/data/exome/processed_final/vcf_processing_all/03_vcf2maf/"



 

process VCF2MAF {
    containerOptions "--bind /data/local/MD_project/data/exome/processed_final/sarek_2_all/variant_calling/mutect2:/data/local/MD_project/data/exome/processed_final/sarek_2_all/variant_calling/mutect2,/data/local/MD_project/data/.vep:/data/local/MD_project/data/.vep,/data/local/MD_project/data/output_vep:/data/local/MD_project/data/output_vep,/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/:/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/,/data/local/MD_project/data/exome/processed_final/vcf_processing_all/02_VEP,/data/local/MD_project/data/exome/processed_final/vcf_processing_all/03_vcf2maf"
    publishDir "${params.outdir}", mode: 'copy'
    maxForks 8

    input:
        path(vcf)

    script:
    """
    #Extract vcf-tumor-id and vcf-normal-id
    vcf_tumor_id=\$(grep -F "##tumor_sample" ${vcf} | cut -d '=' -f2)
    vcf_normal_id=\$(grep -F "##normal_sample" ${vcf} | cut -d '=' -f2)

    /vcf2maf-1.6.22/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf ${vcf} \
        --vcf-tumor-id "\$vcf_tumor_id" \
        --vcf-normal-id "\$vcf_normal_id" \
        --output-maf "${params.outdir}"\$(basename "${vcf}" ".vcf").maf \
        --ref-fasta /data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
        --ncbi-build GRCh38 \
   
    """
}

workflow{
    
    //VEP_ch.view()
    VCF_ch=Channel.fromPath(params.vep_files)
    VCF_ch.view()
    VCF2MAF(VCF_ch)
}