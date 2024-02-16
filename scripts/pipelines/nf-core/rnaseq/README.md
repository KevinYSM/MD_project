# nf-core rnaseq pipeline
<!--
 A bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.
-->

<!--
Sourced from nf-core
 -->

## How to use:
<!--
Just run ./run_command.sh.
If this does not work, copy and paste the command in run_command.sh into the terminal.
You will have to have nextflow installed, and ensure the --input .csv file is correct. This can be generated using the "generate_csv.py" script under nf-core/rnaseq/helper_scripts. 

NB: your PWD (present working directory) has to be the nf-core/rnaseq folder.
-->


### Repository Structure:
<!--
helper_scripts
I have created a helper script called "generate_csv.py".
This script takes a directory with fastqs and generates a .csv file formatted to be run with the nfcore-rnaseq pipeline.
-->

### Test:
<!--
I have not created any tests as of now.
-->


