params.EXOME_raw_reads_directory="/home/ubuntu/scratch/MD_project/data/exome/raw/*/fastqs/*_R{1,2}_*.fastq.gz"



params.trimmed_dir="/home/ubuntu/scratch/MD_project/data/exome/processed/01_trimmed_umis"
params.disambiguated_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/02_disambiguated"


//params.trimmed_umis="/data/local/MD_scholarly/data/processed/exome/01_trimmed_umis/*_R{1,2}_*.fastq.gz"
params.mouse_fasta="/home/ubuntu/scratch/MD_project/data/reference/mouse/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna"
params.human_fasta="/home/ubuntu/scratch/MD_project/data/reference/human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

params.test_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/test"
process AGENT_trim_umis{ 
        maxForks 3
        publishDir params.trimmed_dir, mode: 'copy'
    input:
        file EXOME_read_pair
    output:
        file "*"
        
    """
    java -jar /home/ubuntu/bin/AGeNT/agent/lib/trimmer-3.0.5.jar -fq1 ${EXOME_read_pair[1]} -fq2 ${EXOME_read_pair[2]} -out_loc . -v2
    """
}



workflow{
    read_pairs_ch = Channel.fromFilePairs(params.EXOME_raw_reads_directory, flat: true)
    read_pairs_ch.view()
    //AGENT_trim_umis(read_pairs_ch)


}