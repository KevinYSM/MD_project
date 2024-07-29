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
    --max_memory '120.GB' \
    --input /home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/tumi_tumour.csv \
    --outdir $2 \
    --genome GRCm38 \
    --igenomes_base /data/local/reference/aws/igenomes \
    --max_cpus 63 \
    --wes \
    --save_mapped \
    --skip_tools baserecalibrator,markduplicates


