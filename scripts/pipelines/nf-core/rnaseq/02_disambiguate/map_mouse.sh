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
    --input /home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/trimmed_umis_28122024.csv \
    --outdir $2 \
    --genome GRCm38 \
    --wes \
    --save_mapped \
    --skip_tools baserecalibrator,markduplicates


