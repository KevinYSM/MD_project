params.EXOME_raw_reads_directory="/data/local/MD_scholarly/data/raw/exome/all/renamed/*_R{1,2}_*.fastq.gz"



params.trimmed_dir="/home/ubuntu/scratch/MD_project/data/exome/processed/01_trimmed_umis"
params.disambiguated_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/02_disambiguated"


//params.trimmed_umis="/data/local/MD_scholarly/data/processed/exome/01_trimmed_umis/*_R{1,2}_*.fastq.gz"
params.mouse_fasta="/home/ubuntu/scratch/MD_project/data/reference/mouse/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna"
params.human_fasta="/home/ubuntu/scratch/MD_project/data/reference/human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

params.test_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/test"


process bbsplit{
        maxForks 1
        publishDir params.disambiguated_dir, mode: 'copy'
    input:
        file TRIMMED_read_pair
    output:
        file '*'
    """
    bbsplit.sh in=${TRIMMED_read_pair[1]} \\
     in2=${TRIMMED_read_pair[2]} ref=${params.human_fasta},${params.mouse_fasta} basename=o%.fastq.gz\\
    scafstats=${TRIMMED_read_pair[1]}.scafstats.txt refstats=${TRIMMED_read_pair[1]}.refstats.txt \\
    outu1=${TRIMMED_read_pair[1]}_clean1.fastq.gz outu2=${TRIMMED_read_pair[1]}_clean2.fastq.gz
    """
}

process save_results{
        
        input:
                file FILTERED_fastq
        output:
        """
        echo ${FILTERED_fastq} >> /home/ubuntu/data/local/MD_scholarly/scripts/pipelines/nf-core/sarek/outdir2/a.txt
        cp -s ${FILTERED_fastq[0]} ${FILTERED_fastq[1]} ${FILTERED_fastq[2]} /home/ubuntu/data/local/MD_scholarly/scripts/pipelines/nf-core/sarek/outdir2/
        """
        
}

workflow{
    trimmed_umis_ch=Channel.fromFilePairs(params.trimmed_umis,flat:true)
    bbsplit(trimmed_umis_ch)
    

}