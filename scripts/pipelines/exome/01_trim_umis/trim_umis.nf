params.EXOME_raw_reads_directory="/data/local/MD_scholarly/data/raw/exome/all/renamed/PCB126*_R{1,2}_*.fastq.gz"



params.trimmed_dir="/data/local/MD_scholarly/data/processed/exome/01_trimmed_umis"
params.disambiguated_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/02_disambiguated"


params.trimmed_umis="/data/local/MD_scholarly/data/processed/exome/01_trimmed_umis/*_R{1,2}_*.fastq.gz"
params.mouse_fasta="/data/local/reference/aws/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
params.human_fasta="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

params.test_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/test"
process AGENT_trim_umis{ 
        maxForks 10
        publishDir params.trimmed_dir, mode: 'copy'
    input:
        file EXOME_read_pair
    output:
        file "*"
        
    """
    java -jar $TRIMMER -fq1 ${EXOME_read_pair[1]} -fq2 ${EXOME_read_pair[2]} -out_loc . -v2
    """
}

process bbsplit{
        maxForks 1
        publishDir params.disambiguated_dir, mode: 'copy'
    input:
        file TRIMMED_read_pair
    output:
        file '*'
    """
    /data/local/MD_scholarly/scripts/pipelines/nf-core/sarek/bin/bbmap/bbsplit.sh in=${TRIMMED_read_pair[1]} \\
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
    read_pairs_ch = Channel.fromFilePairs(params.EXOME_raw_reads_directory, flat: true)
    //read_pairs_ch.view()
    trimmed_umis_ch=AGENT_trim_umis(read_pairs_ch).collect()
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