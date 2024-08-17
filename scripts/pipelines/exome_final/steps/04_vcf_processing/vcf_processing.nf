params.vcf_files=params.mutect2
params.outdir="/data/local/MD_project/data/exome/processed_final/vcf_processing"

process VEP {
    containerOptions "--bind /data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2:/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2"
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val(vcf_gz)
    output:
        path("output_vep_updated/*.ann.vcf")

    script:
    """
    vcf=\$(basename "${vcf_gz}" ".vcf.gz").vcf
    bgzip -d -c -f "${vcf_gz}" > "\$vcf"

 

    
    /ensembl-vep/vep \
        --species homo_sapiens \
        --assembly GRCh38 \
        --offline \
        --cache \
        --dir /.vep \
        --input_file /\$(basename "\$vcf") \
        --output_file /output_vep_updated/\$(basename "\$vcf" ".vcf").ann.vcf \
        --everything \
        --vcf \
        --fasta /.vep/homo_sapiens/105_GRCh38/Homo_sapiens_assembly38.fasta \
        --force_overwrite
    """
}

 

process VCF2MAF {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val(vcf)
    output:
        path("*.maf")
    script:
    """
    tumor_id=\$(basename "${vcf}" ".ann.vcf" |  awk -F"_L004" '{print \$1}')
    vcf2maf.pl \
        --inhibit-vep \
        --input-vcf "${vcf}" \
        --output-maf \$(basename "${vcf}" ".vcf").maf \
        --tumor-id "\$tumor_id" \
        --ref-fasta /data/vep_cache/homo_sapiens/105_GRCh38/Homo_sapiens_assembly38.fasta \
        --ncbi-build GRCh38
    """
}

workflow{
    VEP_ch=Channel.fromPath(params.vcf_files)
    VCF_ch=VEP(VEP_ch)
    VCF_ch.view()
    VCF2MAF(VCF_ch)
}