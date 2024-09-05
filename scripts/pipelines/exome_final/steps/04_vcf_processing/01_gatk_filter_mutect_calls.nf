/**Extract from Joakim paper MUM
MuTect 2 calls were
#classified according to quality with GATK FilterMutectCalls and variants failing
those filters were removed. Exceptions were made for known hotspot mutation
sites in GNA11, GNAQ, SF3B1, PLCB4, CYSLTR2, and EIF1AX. Low-quality variants in BAP1 were inspected further for support on DNA and RNA alignments.
Variant annotation was performed with VEP (v. 91.3) and ANNOVAR56 (v. 2016-
05-11), using the databases COSMIC (v. 79), ESP6500 (“siv2_all”), 1000 Genomes
(“2015aug_all”), and dbSNP (“snp138NonFlagged”). For two of the tumors,
exome-sequenced normals were used for further filtering using GATK SelectVariants. Comparisons with TCGA UM mutations were made against mutation lists
downloaded from GDC (accessed May 29, 2017).
*/

params.raw_vcfs="/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2_test/*/*filtered.vcf.gz"
params.outdir="/data/local/MD_project/data/exome/processed_final/vcf_processing/01_filtered"

params.VCF_stats_pair="/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2_test/*/115*sample.mutect2.{filtered.vcf.gz,vcf.gz.stats}"
params.VCF_files = "/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2/*/*.filtered.vcf.gz"
params.VCF_stats="/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2/*/*.vcf.gz.stats"
params.VCF_index="/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2/*/*filtered.vcf.gz.tbi"
process gatk_filter {
    containerOptions "--bind /data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2:/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2,/data/local/MD_project/data/.vep:/data/local/MD_project/data/.vep,/data/local/MD_project/data/output_vep:/data/local/MD_project/data/output_vep,/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta:/data/local/reference/aws/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta,/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai:/data/local/reference/aws/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai,/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict:/data/local/reference/aws/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict"

    input:
        path VCF_trio
    output:
        path "*"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    
    """
    vcf=\$(basename "${VCF_trio[0]}" ".vcf.gz")
    
  
    vcf_basename=\$(basename "${VCF_trio[0]}" ".mutect2.filtered.vcf.gz")
    vcf_dir=/data/local/MD_project/data/exome/processed_final/sarek/variant_calling/mutect2/"\$vcf_basename"/
    tbi_file="\$vcf_dir"/"\$vcf_basename".mutect2.vcf.gz.tbi
    stats_file="\$vcf_dir"/"\$vcf_basename".mutect2.vcf.gz.stats

    cp "${VCF_trio[1]}" "\$vcf_basename".mutect2.filtered.vcf.gz.stats

    echo ${VCF_trio}>a.txt
    /gatk/gatk-4.6.0.0/gatk FilterMutectCalls -V ${VCF_trio[0]} -O  "\$vcf".gatkfiltered.vcf.gz -R "/data/local/reference/aws/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
    --filtering-stats  "${VCF_trio[1]}" --read-index "${VCF_trio[2]}"

    """
}


workflow {
    raw_vcf = Channel.fromPath(params.VCF_files).map { file -> 
        def sample_name = file.baseName
        def replaced_name = sample_name.replaceAll(/\.mutect2.*$/, "")
        
        return [replaced_name, file]
    }
    //raw_vcf.view()

    stats_vcf = Channel.fromPath(params.VCF_stats).map { file -> 
        def sample_name = file.baseName
      
         def replaced_name = sample_name.replaceAll(/\.mutect2.*$/, "")
        
        return [replaced_name, file]
    }

    index_vcf = Channel.fromPath(params.VCF_index).map { file -> 
        def sample_name = file.baseName
      
        def replaced_name = sample_name.replaceAll(/\.mutect2.*$/, "")
        
        return [replaced_name, file]
    }
    //stats_vcf.view()

    vcf_pairs = raw_vcf.join(stats_vcf)
    //vcf_pairs.view()
    //index_vcf.view()
    vcf_trio = vcf_pairs.join(index_vcf).map { it[1..-1] }

    vcf_trio.view()  // To check the files are properly captured by the channel
    //test=Channel.fromPath(params.VCF_files)
    //test.view()
    gatk_filter(vcf_trio)
}