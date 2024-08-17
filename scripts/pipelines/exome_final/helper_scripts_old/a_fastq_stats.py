import os
import re
bbsplit_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/03_bbsplit"
normal_dir="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/02_renamed"


def main():
        samples=get_samples()
        return 0
def get_samples():
        tumour_samples= [file for file in os.listdir(bbsplit_dir) if re.search(".fastq.gz",file)]
        normal_samples= [file for file in os.listdir(normal_dir) if re.search("normal",file)]
        all_samples=tumour_samples+normal_samples
        print(len(tumour_samples))
        print(len(normal_samples))
        print(tumour_samples)
        print(len(all_samples))
 
        return 0
def get_sample_stats():
        return 0
main()