params.renamed_tumis="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/02_renamed/P16*_R{1,2}_*.fastq.gz"
params.mouse_fasta="/data/local/reference/aws/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
params.human_fasta="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

params.bbsplit_dir="/data/local/MD_scholarly/data/processed/exome/03_bbsplit"

process bbsplit{
        maxForks 2
        publishDir params.bbsplit_dir, mode: 'copy'
    input:
        file TRIMMED_read_pair
    output:
        file '*'
    """

    base=\$(echo ${TRIMMED_read_pair[1]} | rev | cut -c 1-9 | rev)
    
    /data/local/MD_scholarly/scripts/pipelines/nf-core/sarek/bin/bbmap/bbsplit.sh in=${TRIMMED_read_pair[1]} \\
     in2=${TRIMMED_read_pair[2]} ref=${params.human_fasta},${params.mouse_fasta} basename=o%.fastq.gz \\
    scafstats=${TRIMMED_read_pair[1]}.scafstats.txt refstats=${TRIMMED_read_pair[1]}.refstats.txt \\
    outu1=${TRIMMED_read_pair[1]}_clean1.fastq.gz outu2=${TRIMMED_read_pair[1]}_clean2.fastq.gz
    """
}

workflow{
    read_pairs_ch = Channel.fromFilePairs(params.renamed_tumis, flat: true)
    //read_pairs_ch.view()
    //tumour_read_pairs_ch=read_pairs_ch.filter(~/.*tumour.*/)
    //bbsplit(tumour_read_pairs_ch)
    bbsplit(read_pairs_ch)
    //trimmed_umis_ch=AGENT_trim_umis(read_pairs_ch).collect()
    //trimmed_umis_ch.view()
    //tumour_read_pairs_ch=trimmed_umis_ch.filter(~/.*tumour.*/)
    //tumour_read_pairs_ch.view()
    //bbsplit(tumour_read_pairs_ch)


    //trimmed_umis_ch=AGENT_trim_umis(read_pairs_ch)
    //need to put a pause here
    
    //trimmed_umis_ch=Channel.fromFilePairs(params.trimmed_umis,flat:true)
    //trimmed_umis_ch.view()
    //output_ch=bbsplit(trimmed_umis_ch)

}