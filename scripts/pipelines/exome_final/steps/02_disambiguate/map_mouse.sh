#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
log_file=logs/$(basename "$0").$(date +"%r%d%h%y" )."$RANDOM".log
exec &> >(tee -a "$log_file")

echo "$CONDA_PREFIX"

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

nextflow run nf-core/sarek \
    -with-singularity $SIF \
    -w $1 \
    -resume \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --input /home/kevin/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/trimmed_samplesheet.csv \
    --outdir $2 \
    --fasta /media/cph/Store4-USB/kevin/references/mouse/ensembl/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --gtf /media/cph/Store4-USB/kevin/references/mouse/ensembl/Mus_musculus.GRCm39.115.gtf \
    --save_reference \
    --wes \
    --save_mapped \
    --skip_tools baserecalibrator,markduplicates


