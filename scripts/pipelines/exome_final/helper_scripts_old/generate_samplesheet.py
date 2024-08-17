import sys
import os
import csv
from os.path import isfile, join
from os import listdir

disambiguated_dir="/data/local/MD_scholarly/data/processed/exome/02_disambiguated"

def get_key( var1, var2):
        num_letters=len(var1)
        counter=0
        for i in range(num_letters):
                if (var1[i]!=var2[i]):
                        break
                counter=counter+1
        return var1[:counter]

def get_fastqs(disambiguated_dir):
    files=[f for f in listdir(disambiguated_dir) if isfile(join(disambiguated_dir,f))]
    fastqs=[]
    print(len(fastqs))

    #Remove all non .fastq.gz files
    for f in files:
        if not (f.endswith(".txt") or f.endswith(".gzscafstats") or f.endswith(".refstats.txt")):
            fastqs.append(f)
    return fastqs




def get_fastq_pairs(fastqs):
    #search given directory
    #take first file
    #use first letter as key, check if only two files match this key
    #if not, use first 2 letters as key, check again
    #repeat adding characters to key until only 2 files match
    #use this number of characters as sample name in output
    fastqs_copy=fastqs
    num_fastqs=len(fastqs)
    corresponding_matches={}
    #for the number of pairs
    for i in range(num_fastqs//2):
        #determine key legnth
        current_fastq=fastqs[0]
        key_string=""
        
        #Determine key string and thus its length
       
        current_key=current_fastq[0:current_fastq.find("_")]

        corresponding_matches[current_key]=[s for s in fastqs if current_key in s]
        print(corresponding_matches)
        for j in corresponding_matches[current_key]:
            fastqs.remove(j) 
    print(num_fastqs)
    print(len(corresponding_matches))
    return corresponding_matches


def createSamplesheet(name):
    if (os.path.exists(name)):
        f=open(name,"w")
        f.truncate()
        f.close()
        return True
    return False


def updateSamplesheet(fastq_pairs, samplesheet, disambiguated_dir):

    f=open(samplesheet,"a")
    writer=csv.writer(f, delimiter="|") 
    writer.writerow(["patient,sex,status,sample,lane,fastq_1,fastq_2"])

    keys_list=fastq_pairs.keys()
    print(len(keys_list))
    i=0
    for sample in keys_list:
        print(sample)
        #assume read 1 is forward
        fastq_pair=fastq_pairs[sample]
        if (("clean1.fastq.gz" in fastq_pair[0])):
            read1=fastq_pair[0]
            read2=fastq_pair[1]
        else:
            read1=fastq_pair[1]
            read2=fastq_pair[0]
        patient=sample[0:sample.find("-")]  
        lane="L001"    
        sex=""  
        status="" 
        row=str(patient+","+sex+","+status+","+sample+","+lane+","+disambiguated_dir+"/"+read1+","+disambiguated_dir+"/"+read2)
        i=i+1
        writer.writerow([row])
        print(row)
        #writer.writerow([str(element+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[0]+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[1]+",auto")])

    print(i)

def main (disambiguated_dir):
    fastqs=get_fastqs(disambiguated_dir)
    paired_fastqs=get_fastq_pairs(fastqs)
    print(paired_fastqs)
    createSamplesheet("/data/local/MD_scholarly/scripts/pipelines/nf-core/sarek/samplesheet.csv")
    updateSamplesheet(paired_fastqs, "/data/local/MD_scholarly/scripts/pipelines/nf-core/sarek/samplesheet.csv",disambiguated_dir)
main(disambiguated_dir)
