params.EXOME_raw_reads_directory="/data/local/MD_project/data/exome/raw/*/fastqs/*_R{1,2}_*.fastq.gz"



params.trimmed_dir="/home/ubuntu/scratch/MD_project/data/exome/processed_final/01_trimmed_umis"


params.test_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/test"
process AGENT_trim_umis{ 
        maxForks 3
        publishDir params.trimmed_dir, mode: 'copy'
    input:
        file EXOME_read_pair
    output:
        file "*"
        
    """
    java -jar /AGeNT/agent/lib/trimmer-3.1.2.jar -fq1 ${EXOME_read_pair[1]} -fq2 ${EXOME_read_pair[2]} -out_loc . -v2
    """
}



workflow{
    read_pairs_ch = Channel.fromFilePairs(params.EXOME_raw_reads_directory, flat: true)
    read_pairs_ch.view()
    AGENT_trim_umis(read_pairs_ch)


}