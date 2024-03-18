Currently having issues with controlfreec as part of Nextflow's Sarek pipeline. Specifically:
"Error: zero reads in windows with the GC-content around 0.35 with interval 0.01, will try again with 0.04
  Number of EM iterations :0
  Error in linear regression, code: -1
  Error in EM => unable to calculate normalized profile
  ERROR: there was a problem in the initial guess of the polynomial. Please contact the support team of change your i
nput parameters. Exit."

Suspect it has something to do with running bbsplit on only tumour samples. This script will output basic statistics for both tumour (disambiguated) and normal samples.