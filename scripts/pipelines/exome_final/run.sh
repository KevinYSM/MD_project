#!/bin/bash
JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/21.0.3-tem"
JAVA_CMD="/home/ubuntu/.sdkman/candidates/java/21.0.3-tem"


nextflow run /home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/steps/01_trim_umis/trim_umis.nf -with-singularity /home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome.sif