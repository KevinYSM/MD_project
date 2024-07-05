#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
log_file=logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}')."$RANDOM".log
exec &> >(tee -a "$log_file")

echo "$CONDA_PREFIX"

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

nextflow run nf-core/sarek \
    -profile singularity \
    -c custom.config \
    -w work_mouse \
    -resume \
    --cf_contamination_adjustment FALSE \
    --cf_contamination 0 \
    --max_memory '120.GB' \
    --input inputs/tumi_tumour.csv \
    --outdir results_mouse \
    --genome GRCm38 \
    --igenomes_base /data/local/reference/aws/igenomes \
    --max_cpus 63 \
    --wes \
    --save_mapped \
    --skip_tools baserecalibrator,markduplicates


