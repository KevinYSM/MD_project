#!/bin/bash
#Build Singularity and bind relevant files
sudo -E env "PATH=$PATH" singularity build /home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome.sif /home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome.def

singularity exec --bind /home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final:/mnt/pipeline
singularity exec --bind /home/ubuntu/data/local/MD_project/data:/mnt/data


singularity run /home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/steps/00_prep/singularity/exome.sif


