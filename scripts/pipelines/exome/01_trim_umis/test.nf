//params.EXOME_raw_reads_directory="/home/ubuntu/scratch/MD_project/data/exome/raw/unbatched/*_R{1,2}_*.fastq.gz"


params.EXOME_raw_reads_directory="data/exome/raw/unbatched/*_R{1,2}_*.fastq.gz"
//params.trimmed_dir="/home/ubuntu/scratch/MD_project/data/exome/processed/01_trimmed_umis"
//params.disambiguated_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/02_disambiguated"



//params.mouse_fasta="/home/ubuntu/scratch/MD_project/data/reference/mouse/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna"
//params.human_fasta="/home/ubuntu/scratch/MD_project/data/reference/human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"


process test_dir{ 
        
    """
    #!/bin/bash
    pwd=\$(pwd | awk -F '/scripts/pipelines/exome/01_trim_umis' '{print \$1}')
    touch \$pwd/${params.EXOME_raw_reads_directory}/apple.txt
    
    """
}



workflow{
    test_dir()
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