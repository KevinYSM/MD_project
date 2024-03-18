cd /home/ubuntu/data/local/MD_scholarly/data/raw/exome || exit
for f in */ ; do
    cd "/home/ubuntu/data/local/MD_scholarly/data/raw/exome/$f" || continue

    cd fastqs/ || continue
    for g in *.fastq.gz; do
        count=$(gunzip -c "$g" | wc -l)
        result="File: $g   Lines: $count"
        echo "$result" >> /home/ubuntu/data/local/MD_scholarly/scripts/pipelines/exome/helper_scripts/output.txt &
    done
    wait  # Wait for all background processes to finish
done