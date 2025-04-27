#!/bin/bash

# Loop through all BAM files in the current directory
for bam in output*.bam; do
    echo "Processing $bam..."

    # Extract header
    samtools view -H "$bam" > header.sam

    # Modify the PL tag in @RG lines
    sed -i 's/PL:[^ \t]*/PL:Illumina/' header.sam

    # Apply the new header
    samtools reheader header.sam "$bam" > "${bam%.bam}.plfixed.bam"

    # Optional: clean up
    rm header.sam

    echo "Modified BAM saved as ${bam%.bam}.plfixed.bam"
done