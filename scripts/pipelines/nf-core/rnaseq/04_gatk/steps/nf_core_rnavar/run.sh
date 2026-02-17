nextflow run nf-core/rnavar -profile singularity --input samplesheet.csv --outdir "/home/user_oruko/work/processed/rna/rnavar_nfcore" --genome GRCh38 \
--dbsnp '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf' \
--dbsnp_tbi '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx' \
--known_indels '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' \
--known_indels_tbi '/home/user_oruko/data/references/Homo_sapiens/GATK_resource_bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'\