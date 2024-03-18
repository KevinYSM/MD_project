import pandas as pd
import os
import re
import shutil
import csv

rna_raw_dir="/home/ubuntu/data/local/MD_scholarly/data/raw/rnaseq/fastqs"

def main():
        fastqs=os.listdir(rna_raw_dir)
        print(len(fastqs))
        fastq_pairs=get_fastq_pairs(fastqs)
      
        
        print(len(fastq_pairs))
        return 0
def get_fastq_pairs(fastqs):
    #search given directory
    #take first file
    #use first letter as key, check if only two files match this key
    #if not, use first 2 letters as key, check again
    #repeat adding characters to key until only 2 files match
    #use this number of characters as sample name in output
    

    files=fastqs
    files_copy=files
    num_fastqs=len(files)
    corresponding_matches={}
    #for the number of pairs
    for i in range(num_fastqs//2):
        #determine key legnth
        current_fastq=files[0]
        key_string=""

        #Determine key string and thus its length
        for j in range(len(current_fastq)):
            key_string=key_string+current_fastq[j]
            matching_values=[s for s in fastqs if key_string in s]
            if (len(matching_values)==2):
                break
        key_length=len(key_string)
        current_key=current_fastq[0:key_length]

     
        corresponding_matches[current_key]=[s for s in fastqs if current_key in s]
        for j in corresponding_matches[current_key]:
            files.remove(j) 
    print(num_fastqs)
    print(len(corresponding_matches))
    return corresponding_matches
def generate_samplesheet():
        return 0

main()