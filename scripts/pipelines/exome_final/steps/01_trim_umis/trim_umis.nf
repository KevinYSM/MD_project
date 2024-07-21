params.EXOME_raw_reads_directory_old="/data/local/MD_project/data/exome/raw/*/fastqs/*_R{1,2}_*.fastq.gz"

params.EXOME_raw_reads=  "/data/local/MD_project/data/exome/raw/*/fastqs/PCB-31-TIL*_R{1,2}_*.fastq.gz"

params.UNFINISHED="/data/local/MD_project/data/exome/raw/*/PCB-31-TIL*.fastq.gz"



process AGENT_trim_umis{ 
        maxForks 3
        publishDir params.TRIMMED_DIR, mode: 'copy'
    input:
        file EXOME_read_pair
    output:
        file "*"
        
    """
    java -jar /AGeNT/agent/lib/trimmer-3.1.2.jar -fq1 ${EXOME_read_pair[1]} -fq2 ${EXOME_read_pair[2]} -out_loc . -v2
    """
}



workflow{
    read_pairs_ch = Channel.fromFilePairs(params.EXOME_raw_reads, flat: true)
    read_pairs_ch.view()
    AGENT_trim_umis(read_pairs_ch)
}