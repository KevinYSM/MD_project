#!/bin/bash
SIF="/home/user_oruko/work/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome_2.sif"
nextflow run bbsplit.nf \
-c exome.config \
/
