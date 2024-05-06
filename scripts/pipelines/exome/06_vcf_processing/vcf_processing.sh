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

JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/current"
JAVA_CMD="/home/ubuntu/.sdkman/candidates/java/current/bin/java"
SINGULARITY_TMPDIR="/data/local/MD_project/scripts/pipelines/exome/05_sarek_joakim/pipeline_tmp/s_tmp"
SINGULARITY_LOCALCACHEDIR="/data/local/MD_project/scripts/pipelines/exome/05_sarek_joakim/pipeline_tmp/s_cache"

nextflow run vcf_processing.nf \
    -resume \
    --max_memory '120.GB' \
    --max_cpus 63 