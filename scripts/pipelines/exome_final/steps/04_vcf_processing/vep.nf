params.vcf_files="/home/user_oruko/work/processed/exome/sarek_intervals/variant_calling/mutect2/*/*filtered.vcf.gz"
params.outdir="/home/user_oruko/work/processed/exome/vcf_processing_intervals"

process VEP {
    containerOptions "--bind /home/user_oruko/work/processed/exome/sarek_intervals/variant_calling/mutect2:/home/user_oruko/work/processed/exome/sarek_intervals/variant_calling/mutect2,/home/user_oruko/work/processed/exome/.vep:/home/user_oruko/work/processed/exome/.vep,/home/user_oruko/work/processed/exome/vcf_processing_intervals/output_vep:/home/user_oruko/work/processed/exome/vcf_processing_intervals/output_vep,/home/user_oruko/data/references/aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta:/home/user_oruko/data/references/aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta,/home/user_oruko/work/processed/exome/vcf_processing_intervals/output_vep_refseq/:/home/user_oruko/work/processed/exome/vcf_processing_intervals/output_vep_refseq/"
    publishDir "${params.outdir}", mode: 'copy'
    maxForks 1

    input:
        val(vcf_gz)
    

    script:
    """
    vcf=\$(basename "${vcf_gz}" ".vcf.gz").vcf
    bgzip -d -c -f "${vcf_gz}" > "\$vcf"

 

    
    /ensembl-vep/vep \
        --species homo_sapiens \
        --use_given_ref \
        --assembly GRCh38 \
        --refseq \
        --offline \
        --cache \
        --dir /home/user_oruko/work/processed/exome/.vep \
        --dir_cache /home/user_oruko/work/processed/exome/.vep \
        --input_file "\$vcf" \
        --output_file /home/user_oruko/work/processed/exome/vcf_processing_intervals/output_vep_refseq/\$(basename "\$vcf" ".vcf").ann.vcf \
        --everything \
        --vcf \
        --fork 2 \
        --fasta /home/user_oruko/data/references/aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
        --force_overwrite
    """
}

 



workflow{
    VEP_ch=Channel.fromPath(params.vcf_files)
    //VEP_ch.view()
    VCF_ch=VEP(VEP_ch)
}