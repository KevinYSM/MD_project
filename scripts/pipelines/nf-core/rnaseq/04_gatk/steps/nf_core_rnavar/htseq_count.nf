#!/usr/bin/env nextflow

// Define parameters
params.bam_dir = "/home/user_oruko/work/processed/rna/star_salmon/*.bam"  // Directory containing the BAM files
params.gtf_file = "/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.108.gtf"             // Path to your GTF annotation file.  **REQUIRED**
params.strandedness = "no"       // Strandedness of your library: "yes", "no", or "reverse". **REQUIRED**
params.id_attribute = "gene_id"  // GTF attribute to use for counting (e.g., "gene_id", "transcript_id")
params.out_dir = "/home/user_oruko/work/processed/rna_other/htseq_count"    // Directory for HTSeq-count output


// Process to run HTSeq-count
process htseq_count {

    input:
    path bam

    output:


    
    
    script:
    

    """
    base=\$(basename ${bam})

    htseq-count \
        -f bam \
        -r pos \
        -s no \
        -a 10 \
        -t exon \
        -i gene_id \
        -m intersection-nonempty \
        ${bam} \
        ${params.gtf_file} > ${params.out_dir}/\${base}.htseq.counts


    
    """
    
    
    
}

workflow {

    // Get the list of BAM files
    bam_ch=Channel.fromPath(params.bam_dir)
    htseq_count(bam_ch)
    
}
