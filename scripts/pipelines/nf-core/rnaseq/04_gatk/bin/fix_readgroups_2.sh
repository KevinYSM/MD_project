#!/bin/bash

# === CONFIG ===
INPUT_DIR="/home/user_oruko/work/processed/rna/gatk/BAM_merge_alignment"              # Change to your input BAM directory
OUTPUT_DIR="/home/user_oruko/work/processed/rna/gatk/BAM_replaced_read_groups"       # Change to your desired output directory


# === PROCESS BAM FILES ===
for BAM in "$INPUT_DIR"/*.bam; do
    BASENAME=$(basename "$BAM" .bam)
    FIXED_BAM="${OUTPUT_DIR}/${BASENAME}.fixed.bam"

    patient_id="${BASENAME%%_m*}"  # Extract everything before '_p'
    patient_id="${patient_id#p}"   # Remove the leading 'p'
    echo "Processing $BAM → $FIXED_BAM"

    # Replace incorrect RG:Z:57_SXXX tags with RG:Z:SXXX (match from actual header ID)
    HEADER_RG_ID=$(samtools view -H "$BAM" | grep '^@RG' | sed -n 's/.*ID:\([^\t]*\).*/\1/p')
    

    HEADER_RG=$(samtools view -H "$BAM" | grep '^@RG')


    declare -a patient_ids_152=("9" "16" "20" "28" "30" "31" "36" "39" "50" "54")

    
   
    sample_group="S301"    
    for id in ${patient_ids_152[@]}; do
        
        if [[ ${patient_id} == "${id}" ]]; then
            sample_group="S152"
            
        fi
    done
    

    NEW_HEADER="@RG\tID:${patient_id}_${sample_group}\tSM:${patient_id}\tLB:${sample_group}\tPL:Illumina"
    # Write the header to a temporary file
    echo -e "$NEW_HEADER" > "$OUTPUT_DIR/header.sam"
    # Step 1: Reheader the BAM file and save it to a temporary file
    samtools reheader <(echo -e "$NEW_HEADER") "$BAM" | \
        samtools view -h | \
        sed "s/RG:Z:${sample_group}/RG:Z:${patient_id}_${sample_group}/g" | \
        samtools view -bS -o "$FIXED_BAM" 

    echo "Indexing $FIXED_BAM"
    samtools index "$FIXED_BAM"

    echo "✅ Done: $FIXED_BAM"
done

echo "🎉 All BAM files processed and saved in: $OUTPUT_DIR"

