#!/bin/bash

# === CONFIG ===
INPUT_DIR="/home/user_oruko/work/processed/rna/gatk/BAM_merge_alignment"              # Change to your input BAM directory
OUTPUT_DIR="/home/user_oruko/work/processed/rna/gatk/BAM_replaced_read_groups_2"       # Change to your desired output directory


# === PROCESS BAM FILES ===
for BAM in "$INPUT_DIR"/*.bam; do
    BASENAME=$(basename "$BAM" .bam)
    FIXED_BAM="${OUTPUT_DIR}/${BASENAME}.fixed.bam"
    FINAL_BAM="${OUTPUT_DIR}/${BASENAME}.final.bam"

    patient_id="${BASENAME%%_m*}"  # Extract everything before '_p'
    patient_id="${patient_id#p}"   # Remove the leading 'p'
    echo "Processing $BAM → $FIXED_BAM"

    # Replace incorrect RG:Z:57_SXXX tags with RG:Z:SXXX (match from actual header ID)
    HEADER_RG_ID=$(samtools view -H "$BAM" | grep '^@RG' | sed -n 's/.*ID:\([^\t]*\).*/\1/p')
    

    

    declare -a patient_ids_152=("9" "16" "20" "28" "30" "31" "36" "39" "50" "54")

    
   
    sample_group="S301"    
    for id in ${patient_ids_152[@]}; do
        
        if [[ ${patient_id} == "${id}" ]]; then
            sample_group="S152"
            
        fi
    done
    
    # Update the RG tags in the reads: Replace incorrect RG:Z tags with the correct ID
    (
    NEW_HEADER="@RG\tID:${patient_id}_${sample_group}\tSM:${patient_id}\tLB:${sample_group}\tPL:Illumina"
    #samtools reheader 

    samtools view -h "$BAM" | \
        sed -e "s/RG:Z:${sample_group}/RG:Z:${patient_id}_${sample_group}/g" -e "s/@RG\tID:${sample_group}\tSM:${patient_id}_${samplegroup}/${NEW_HEADER}/g"| \
        samtools view -bS -o "$FIXED_BAM" - 

    

    echo "Indexing $FINAL_BAM"
    #samtools index "$FINAL_BAM"
    ) &
    wait
    echo "✅ Done: $FINAL_BAM"
    
done

wait

echo "🎉 All BAM files processed and saved in: $OUTPUT_DIR"

