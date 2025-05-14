//params.EXOME_raw_reads_directory_old="/data/local/MD_project/data/exome/raw/*/fastqs/*_R{1,2}_*.fastq.gz"

params.RNA_raw_reads=  "/home/user_oruko/data/raw/rna/all_fastqs/*_R{1,2}_*.fastq.gz"



process AGENT_trim_umis{ 
        maxForks 5
        publishDir params.TRIMMED_DIR, mode: 'copy'
    input:
        file RNA_read_pair
    output:
        file "*"
        
    """
    /AGeNT/agent/agent.sh trim -fq1 ${RNA_read_pair[1]} -fq2 ${RNA_read_pair[2]} -v2 -out_loc .
    """
}



workflow{
    read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads, flat: true)
    read_pairs_ch.view()
    AGENT_trim_umis(read_pairs_ch)
}