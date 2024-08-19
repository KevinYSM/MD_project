#!/bin/bash

JAVA_HOME="/home/ubuntu/.sdkman/candidates/java/21.0.2-tem"
JAVA_CMD="/home/ubuntu/.sdkman/candidates/java/21.0.2-tem"
SINGULARITY_TMPDIR="/data/local/MD_project/scripts/pipelines/exome/05_sarek_joakim/pipeline_tmp/s_tmp"
SINGULARITY_LOCALCACHEDIR="/data/local/MD_project/scripts/pipelines/exome/05_sarek_joakim/pipeline_tmp/s_cache"
source /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/02_disambiguate_joakim/map_human.sh
source /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/02_disambiguate_joakim/map_mouse.sh