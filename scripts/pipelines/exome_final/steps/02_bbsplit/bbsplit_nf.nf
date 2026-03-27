#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include the nf-core module logic
include { BBMAP_BBSPLIT } from './modules/nf-core/bbmap/bbsplit/main'

params.TUMI_raw_reads = "/media/cph/Store4-USB/kevin/exome/data/raw/*_R{1,2}_*.fastq.gz"
params.bbsplit_index   = "/media/cph/Store4-USB/kevin/references/bbsplit_human_mouse_index/"

process AGENT_trim_umis { 
    tag "${meta.id}"
    maxForks 5
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("*.fastq.gz"), emit: trimmed_reads
        
    script:
    // AGeNT usually outputs files with specific suffixes; adjusting for your script
    """
    /opt/AGeNT/agent/agent.sh trim \
        -fq1 ${reads[0]} \
        -fq2 ${reads[1]} \
        -v2 -out_loc .
    """
}

workflow {
    // 1. Setup Input Channel
    // nf-core modules expect a 'meta' map for sample tracking
    read_pairs_ch = Channel.fromFilePairs(params.EXOME_raw_reads, checkIfExists: true)
        .map { id, files -> [ [id:id, single_end:false], files ] }

    // 2. Trim UMIs
    AGENT_trim_umis(read_pairs_ch)

    // 3. Disambiguate using BBSplit
    // We pass 'false' for build_index because we are pointing to a pre-built path
    BBMAP_BBSPLIT (
        AGENT_trim_umis.out.trimmed_reads,
        params.bbsplit_index,
        false 
    )

    // 4. Access the Human-only reads for Sarek
    // BBSplit module outputs a list of files; you'll want the one matching the human reference
    BBMAP_BBSPLIT.out.primary_reads.view { meta, fastqs -> "Human reads for ${meta.id}: ${fastqs}" }
}