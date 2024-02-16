#a script to generate the samplesheet required to run nf-core's RNASeq pipeline

#Input: directory with all RAW fastq files
#Output: samplesheet.csv


import sys
import os
import csv



def createSamplesheet(name):
    if (os.path.exists(name)):
        f=open(name,"w")
        f.truncate()
        f.close()
        return True
    return False

def invalidArgs():
    if (len(sys.argv)==1):
        return True
    return False

def updateSamplesheet(add, name):
    fastq_pairs=get_fastq_pairs()
    f=open(name,"a")
    writer=csv.writer(f, delimiter="|") 
    if (add!=True):
        writer.writerow(["sample,fastq_1,fastq_2,strandedness"])

    keys_list=fastq_pairs.keys()
    print(len(keys_list))
    i=0
    for element in keys_list:
        print(element)
        #assume read 1 is forward
        fastq_pair=fastq_pairs[element]
        if not (("R1_001.fastq.gz" in fastq_pair[0]) or ("_1.fastq.gz" in fastq_pair[0])):
            continue
       
        
        row=str(element+","+fastq_pair[0]+","+fastq_pair[1]+",auto")
        i=i+1
        writer.writerow([str(element+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[0]+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[1]+",auto")])

    print(i)
def get_fastq_pairs():
    #search given directory
    #take first file
    #use first letter as key, check if only two files match this key
    #if not, use first 2 letters as key, check again
    #repeat adding characters to key until only 2 files match
    #use this number of characters as sample name in output
    
    if (len(sys.argv)!=4):
        print("Directory or samplesheet name not given!")
        return 0
    fastq_directory=sys.argv[3]
    files=os.listdir(fastq_directory)
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
            matching_values=[s for s in os.listdir(fastq_directory) if key_string in s]
            if (len(matching_values)==2):
                break
        key_length=len(key_string)
        current_key=current_fastq[0:key_length]

     
        corresponding_matches[current_key]=[s for s in os.listdir(fastq_directory) if current_key in s]
        for j in corresponding_matches[current_key]:
            files.remove(j) 
    print(num_fastqs)
    print(len(corresponding_matches))
    return corresponding_matches



    
def main():
   
    invalid_option=False
    if (invalidArgs()):
        print("Invalid option.")
        print("Try 'python3 generate_csv.py -h' for more information")
        return 0

    flag=sys.argv[1]
    
    if (flag=="-h"):
        print("Usage:\tpython3\tgenerate_csv.py\t[FLAG]\t[DIRECTORY]\n")
        print("-g/--generate:\tgenerate a new samplesheet.csv (deleting old one) using given directory fastqs")
        print("-a/--add:\add to existing samplesheet.csv (deleting old one) using given directory fastqs")
       
        print("\nFor use with nfcore RNASeq pipeline")
        return 0
    

    elif (flag=="-g" or flag=="--generate"):
        print("Generating new samplesheet")
        samplesheet_name=sys.argv[2]
        createSamplesheet(samplesheet_name)
        updateSamplesheet(False, samplesheet_name)
        print(samplesheet_name +" was successfully created!")

    elif (flag=="-a" or flag=="--add"):
        samplesheet_name=sys.argv[2]
        print("Adding to samplesheet.csv")
        if not (os.path.exists(samplesheet_name)):
            print("Can't add to a samplesheet that does not exist!")
            return 0
        updateSamplesheet(True, samplesheet_name)
        print(samplesheet_name +" was successfully updated!")
                
   
    else:
        print("Invalid option.")
        print("Try 'python3 generate_csv.py -h' for more information")
    return 0


test_data_directory="/home/ubuntu/MD_scholarly/data/test/GSE110004/"
main()