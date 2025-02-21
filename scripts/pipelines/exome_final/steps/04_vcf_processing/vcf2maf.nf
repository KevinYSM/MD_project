params.vep_files="/home/user_oruko/work/processed/exome/vcf_processing/filtered_vcfs/*filtered.vcf"
params.outdir="/home/user_oruko/work/processed/exome/vcf_processing/output_maf"



 

process VCF2MAF {
    containerOptions "--bind /home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2:/home/user_oruko/work/processed/exome/sarek/variant_calling/mutect2,/home/user_oruko/work/processed/exome/.vep:/home/user_oruko/work/processed/exome/.vep,/home/user_oruko/work/processed/exome/vcf_processing/output_vep:/home/user_oruko/work/processed/exome/vcf_processing/output_vep,/home/user_oruko/data/references/aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta:/home/user_oruko/data/references/aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta,/home/user_oruko/work/processed/exome/vcf_processing/output_vep_refseq/:/home/user_oruko/work/processed/exome/vcf_processing/output_vep_refseq/"
    publishDir "${params.outdir}", mode: 'copy'
    maxForks 8

    input:
        path(vcf)

    script:
    """
    vcf=${vcf}


    #Extract vcf-tumor-id and vcf-normal-id
    vcf_tumor_id=\$(grep -F "##tumor_sample" "\${vcf}" | cut -d '=' -f2)
    vcf_normal_id=\$(grep -F "##normal_sample" "\${vcf}" | cut -d '=' -f2)

    /vcf2maf-1.6.22/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf "\${vcf}" \
        --vcf-tumor-id "\$vcf_tumor_id" \
        --vcf-normal-id "\$vcf_normal_id" \
        --output-maf /home/user_oruko/work/processed/exome/vcf_processing/output_maf/\$(basename "\${vcf}" ".vcf").maf \
        --ref-fasta /home/user_oruko/data/references/aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
        --ncbi-build GRCh38 \
   
    """
}

workflow{
    
    //VEP_ch.view()
    VCF_ch=Channel.fromPath(params.vep_files)
    VCF_ch.view()
    VCF2MAF(VCF_ch)
}