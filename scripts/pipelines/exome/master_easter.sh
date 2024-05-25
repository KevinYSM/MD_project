#!/bin/bash
source /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/02_disambiguate_test/run_command.sh
python3 /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/helper_scripts/04_generate_samplesheets_copy.py
source /home/ubuntu/data/local/MD_project/scripts/pipelines/exome/05_sarek/run_command.sh