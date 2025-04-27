#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam="$1"
#output_dir="$2"
#gtf=genome/Homo_sapiens.GRCh37.75.sorted.gtf
gtf=/home/user_oruko/data/references/Homo_sapiens/i/Homo_sapiens.GRCh38.108.gtf

htseq-count -r name -q -f bam -s reverse -m intersection-strict "$bam" "$gtf" > $(basename "$bam" ".bam").s_reverse.gene_counts