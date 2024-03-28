#!/bin/bash
pwd=$(pwd | awk -F '/scripts/pipelines/exome/00_collate_and_rename' '{print $1}')

python3 $pwd/scripts/pipelines/exome/helper_scripts/00_collate_and_rename.py $pwd/data/exome/raw/batched