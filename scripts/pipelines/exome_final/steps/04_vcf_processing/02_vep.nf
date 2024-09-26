params.vcf_files_old="/data/local/MD_project/data/exome/processed_final/vcf_processing_all/01_filtered/*vcf.gz"
params.vcf_files="/data/local/MD_project/data/exome/processed_final/sarek_2_all/variant_calling/mutect2/*/*.filtered.vcf.gz"
params.outdir="/data/local/MD_project/data/exome/processed_final/vcf_processing_all/02_VEP"

process VEP {
    containerOptions "--bind /data/local/MD_project/data/exome/processed_final/sarek_2_all/variant_calling/mutect2:/data/local/MD_project/data/exome/processed_final/sarek_2_all/variant_calling/mutect2,/data/local/MD_project/data/.vep:/data/local/MD_project/data/.vep,/data/local/MD_project/data/output_vep:/data/local/MD_project/data/output_vep,/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/:/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/,/data/local/MD_project/data/exome/processed_final/vcf_processing_all/01_filtered/,/data/local/MD_project/data/exome/processed_final/vcf_processing_all/02_VEP"
    publishDir "${params.outdir}", mode: 'copy'
    maxForks 4

    input:
        val(vcf_gz)
    

    script:
    """
    vcf=\$(basename "${vcf_gz}" ".vcf.gz").vcf
    bgzip -d -c -f "${vcf_gz}" > "\$vcf"

 

    
    /ensembl-vep/vep \
        --species homo_sapiens \
        --assembly GRCh38 \
        --offline \
        --cache \
        --dir /data/local/MD_project/data/.vep \
        --input_file "\$vcf" \
        --output_file "${params.outdir}"/\$(basename "\$vcf" ".vcf").ann.vcf \
        --everything \
        --vcf \
        --format vcf \
        --fork 4 \
        --fasta /data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
        --force_overwrite
    """
}

 



workflow{
    VEP_ch=Channel.fromPath(params.vcf_files)
    //VEP_ch.view()
    VCF_ch=VEP(VEP_ch)
}