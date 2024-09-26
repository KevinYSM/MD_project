params.trimmed_umis="/data/local/MD_project/data/exome/processed_final/01_trimmed_umis_redone_/*{R1,R2}*.fastq.gz"
params.mouse_fasta=params.FASTA_HUMAN
params.human_fasta=params.FASTA_MOUSE
params.list = '/data/local/MD_project/scripts/pipelines/exome_final/steps/02_bbsplit/tumour_paths.txt'


process test_rename{
    publishDir params.disambiguated_dir, mode: 'copy'
    input:
        file TRIMMED_read_pair
    output:
        file '*'
    """
    rename_1=\$(basename \$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_rename.py ${TRIMMED_read_pair[1]}) .fastq.gz)
    rename_2=\$(basename \$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_rename.py ${TRIMMED_read_pair[2]}) .fastq.gz)


    batch=\$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_batch.py ${TRIMMED_read_pair[1]})

    if  [[ "\${batch}" == "374" ]]
    then
    echo ${TRIMMED_read_pair[1]} >> \${rename_1}.txt
    touch b.txt
    fi
    touch a.txt 
    touch "\${batch}".txt

    
    """
}

process bbsplit_batch{
        maxForks 2
        
        publishDir params.disambiguated_dir, mode: 'copy'
    input:
        file TRIMMED_read_pair
    output:
        file '*'
    """
    rename_1=\$(basename \$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_rename.py ${TRIMMED_read_pair[1]}) .fastq.gz)
    rename_2=\$(basename \$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_rename.py ${TRIMMED_read_pair[2]}) .fastq.gz)
    batch=\$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_batch.py ${TRIMMED_read_pair[1]})


    if  [[ "\${batch}" == "374" ]]
    then
    touch a.txt
    

    fi

    
    echo ${TRIMMED_read_pair[1]} >> \${rename_1}.txt
    

    bbsplit.sh -Xmx50g in=${TRIMMED_read_pair[1]} \\
    in2=${TRIMMED_read_pair[2]} ref=${params.human_fasta},${params.mouse_fasta} basename=o%.fastq.gz\\
    scafstats=\${rename_1}.scafstats.txt refstats=\${rename_1}.refstats.txt \\
    outu1=\${rename_1}_clean1.fastq.gz outu2=\${rename_2}_clean2.fastq.gz

    """
}

process bbsplit{
        maxForks 2
        
        publishDir params.BBSPLIT_DIR, mode: 'copy'
    input:
        file TRIMMED_read_pair
    output:
        file '*'
    """
    

    rename=\$(basename ${TRIMMED_read_pair[0]} .fastq.gz)

    
    bbsplit.sh -Xmx50g in=${TRIMMED_read_pair[0]} in2=${TRIMMED_read_pair[1]} \
    ref_human=${params.human_fasta} ref_mouse=${params.mouse_fasta} \
    basename=out_%_\${rename}_#.fastq.gz scafstats=\${rename}.scafstats.txt refstats=\${rename}.refstats.txt \
    ambiguous2=all \
    path=/data/local/MD_project/scripts/pipelines/exome_final/ref
"""
}

process bbsplit_old{
        maxForks 2
        
        publishDir params.BBSPLIT_DIR, mode: 'copy'
    input:
        file TRIMMED_read_pair
    output:
        
    """
    rename_1=\$(basename \$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_rename.py ${TRIMMED_read_pair[1]}) .fastq.gz)
    rename_2=\$(basename \$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_rename.py ${TRIMMED_read_pair[2]}) .fastq.gz)
    batch=\$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_batch.py ${TRIMMED_read_pair[1]})
    tumour=\$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_tumour.py ${TRIMMED_read_pair[1]})

    rename=\$(basename ${TRIMMED_read_pair[1]} .fastq.gz)

    tumour=\$(python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/return_tumour.py ${TRIMMED_read_pair[1]})
    echo ${TRIMMED_read_pair[1]} \${tumour} >> /data/local/MD_project/scripts/pipelines/exome_final/steps/02_bbsplit/a.txt
    if  [[ "\${tumour}" == "tumour" ]]
    then
        bbsplit.sh -Xmx50g in=${TRIMMED_read_pair[1]} \
        in2=${TRIMMED_read_pair[2]} ref_human=${params.human_fasta} ref_mouse=${params.mouse_fasta} \
        basename=out_%_\${rename_1}_#.fastq.gz scafstats=\${rename_1}.scafstats.txt refstats=\${rename_1}.refstats.txt \
        minid=0.9 ambiguous2=all \
        path=/data/local/MD_project/scripts/pipelines/exome_final/
    
    fi
    
   

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
    trimmed_umis_ch = Channel
    .fromPath(params.list)
    .splitText() { it.trim() }
    .map { file(it) }
    .collect()
    .map { files -> 
        def paired = files.collate(2) // Groups the files into pairs assuming they are listed sequentially
        return paired.collect { pair -> 
            def r1 = pair.find { it.toString().contains('_R1_') }
            def r2 = pair.find { it.toString().contains('_R2_') }
            return [r1, r2]
        }
    } .flatMap { it }
   

    trimmed_umis_ch.view()

    //trimmed_umis_ch=Channel.fromFilePairs(params.trimmed_umis,flat:true)
    
    //test_rename(trimmed_umis_ch)
    bbsplit(trimmed_umis_ch)
    

}