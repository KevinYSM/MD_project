#!/usr/bin/env bash
echo "Starting!"
python3 /data/local/MD_project/scripts/pipelines/exome_final/ngs_disambiguate/disambiguate/disambiguate.py -s "109_tumour_sample" \
-o "./" -a bwa "/data/local/MD_project/scripts/pipelines/exome_final/work/41/16b118bf1114d5652ab6025571250f/109_tumour_sample.human.bam" \
"/data/local/MD_project/scripts/pipelines/exome_final/work/41/16b118bf1114d5652ab6025571250f/109_tumour_sample.mouse.bam" \
