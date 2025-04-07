//Generate unaligned bam files
process UnalignedBams{
        maxForks 5
        publishDir '/home/user_oruko/work/processed/rna/gatk/BAM_unaligned', mode: 'copy'

        
        input:
        tuple val(sample_id), path(read_files)
        

        output:
        path '*.bam'
        script:
        """
        #Patient IDs with read group name S152
        declare -a patient_ids_152=("9" "16" "20" "28" "30" "31" "36" "39" "50" "54")

        
        patient_id=\$(echo "${read_files[0]}" | sed 's/p\\([0-9]\\{1,3\\}\\)_primary.*\$/\\1/')
        sample_group="S301"    
        for id in \${patient_ids_152[@]}; do
            
            if [[ \${patient_id} == "\${id}" ]]; then
                sample_group="S152"
                
            fi
        done


        java -jar /home/user_oruko/other/picard/picard.jar FastqToSam \
        F1=${read_files[0]} \
        F2=${read_files[1]} \
        O=unaligned_"${sample_id}".bam \
        SM=\${patient_id} \
        READ_GROUP_NAME=\${sample_group} \
        PLATFORM=ILLUMINA
        
        """
}
//Data cleanup using MergeBamAlignment and Markduplicates
process MergeBamAlignment{
    input:
        file BAM_file
        

    output:
        path "*.bam"

        """
        base=\$(basename $BAM_file)
        java -jar /picard/picard.jar MergeBamAlignment \
        --ALIGNED_BAM $BAM_file \\
        UNMAPPED=unmapped.bam \\
        O=merge_alignments.bam \\
        R=reference_sequence.fasta

        
        """
    
}