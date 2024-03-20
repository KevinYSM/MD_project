#!/bin/bash
#NB: all raw exome data (.fastq.gz) needs to be in a single folder
raw_exome_data="/home/ubuntu/scratch/MD_project/data/exome/raw/unbatched"
metadata_table=
human_reference=
mouse_reference=
AGeNT_loc=

#00_collate and rename
./home/ubuntu/scratch/MD_project/scripts/pipelines/exome/00_collate_and_rename/run_command.sh  $raw_exome_data
#01_trim_umis
./home/ubuntu/scratch/MD_project/scripts/pipelines/exome/01_trim_umis/run_command.sh

#02_rename
#03_disambiguate
#04_samplesheets
#05_sarek