Takes PDX fastq files and removes mouse reads. 

Important:
1. Do not use BBSplit on TILs (can attempt this later for interests sake)


/data/local/MD_scholarly/scripts/pipelines/nf-core/sarek/bin/bbmap/bbsplit.sh in=${TRIMMED_read_pair[1]} \\
     in2=${TRIMMED_read_pair[2]} ref=${params.human_fasta},${params.mouse_fasta} basename=o%.fastq.gz\\
    scafstats=${TRIMMED_read_pair[1]}.scafstats.txt refstats=${TRIMMED_read_pair[1]}.refstats.txt \\
    outu1=${base}_clean1.fastq.gz outu2=${base}_clean2.fastq.gz