params.vcf_files="/home/ubuntu/scratch/MD_project/data/exome/processed/variant_calling_0605/04_sarek_joakim_disambiguate/annotation/mutect2/*/*_VEP.ann.vcf.gz"
params.outdir="/home/ubuntu/scratch/MD_project/data/exome/processed/maf/"

process VEP {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val(vcf_gz)
    output:
        path("output_vep_updated/*.ann.vcf")

    script:
    """
    vcf=\$(basename "${vcf_gz}" ".vcf.gz").vcf
    bgzip -dcf "${vcf_gz}" > "\$vcf"

    echo \$vcf

    if [ ! -d output_vep_updated ]
    then
        mkdir output_vep_updated

    fi
    vep \
        --species homo_sapiens \
        --assembly GRCh38 \
        --input_file \$vcf \
        --output_file output_vep_updated/\$(basename "\$vcf" ".vcf").ann.vcf \
        --everything \
        --vcf \
        --cache \
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
    perl /home/ubuntu/bin/mskcc-vcf2maf-754d68a/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf "${vcf}" \
        --output-maf \$(basename "${vcf}" ".vcf").maf \
        --tumor-id "\$tumor_id" \
        --ref-fasta /home/ubuntu/scratch/MD_project/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        --ncbi-build GRCh38
    """
}

workflow{
    VEP_ch=Channel.fromPath(params.vcf_files)
    //VCF_ch=VEP(VEP_ch)
    //VCF_ch.view()
    VCF2MAF(VEP_ch)
}