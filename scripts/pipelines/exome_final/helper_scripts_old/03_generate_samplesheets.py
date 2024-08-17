import pandas as pd
import os
import re
import shutil
import csv

tumour_dir="/data/local/MD_project/data/exome/processed_final/02_disambiguated/results_disambiguate/bam_to_fastq"
normal_dir="/data/local/MD_scholarly/data/processed/exome/02_renamed"
batch_ids=["335","374","693","952"]

def main():
        tumour_samples= os.listdir(tumour_dir)
        return 0

def get_batches(filtered_files, batch_ids):
        batches=[]

        for batch in batch_ids:
                filter=re.compile(re.escape("B"+batch))
                batch_files=[file for file in filtered_files if filter.search(file)]
                batches.append(batch_files)
        return batches


def get_fastqs(bbsplit_dir):
        tumour_samples= os.listdir(bbsplit_dir)
        
        normal_samples= [file for file in os.listdir(normal_dir) if re.search("normal",file)]
     
        all_samples=tumour_samples+normal_samples
        filter=re.compile(r'\.fastq\.gz$')
        filtered_files=[file for file in all_samples if filter.search(file)]
        print(len(filtered_files))

        all_csv="/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/all.csv"
        if (os.path.exists(all_csv)):
                f=open(all_csv,"w")
                f.truncate()
                f.close()

        batches=get_batches(filtered_files,batch_ids)
        
        generate_samplesheets(batches)

        sort_samplesheet("/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/all.csv")

        return 0

def generate_samplesheets(batches): 
        paired={}
        create_all_samplesheet(batches)
        #for batch in batches:
        #      create_paired_samplesheet(batch)
              
           


def create_paired_samplesheet(batch):
       
        #filter samplesheet
        paired={}

        for f in batch:
                
                attributes_f=get_attributes(f)
                for g in batch:
                        attributes_g=get_attributes(g)
         
                        if (attributes_f["read"]=="1"):

                                if (attributes_f["ID"]==attributes_g["ID"]
                                    and attributes_f["status"]==attributes_g["status"]
                                    and attributes_g["read"]=="2"):
                                        pair=[f,g]
                                        paired[attributes_f["ID"]]=pair
                                
                                      
               

                                      
        #update_samplesheet(batch,paired_name)
  
        #Create samplesheet   
        batch_id=get_attributes(batch[0])["batch"]                          
        paired_name="/home/ubuntu/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/"+batch_id+"_paired.csv"
        if (os.path.exists(paired_name)):
                f=open(paired_name,"w")
                f.truncate()
                f.close()
        
        #Update samplesheet
        update_samplesheet(paired,paired_name)
        #update_samplesheet(paired,"/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/all.csv")
        return paired

def create_all_samplesheet(batches):
        all_csv="/home/ubuntu/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/"+"all.csv"
        if (os.path.exists(all_csv)):
                f=open(all_csv,"w")
                f.truncate()
                f.close()
        paired=[]
        pair_count=0
        for batch in batches:
                batch_id=get_attributes(batch[0])["ID"]
                
                for f in batch:
                
                        attributes_f=get_attributes(f)
                        for g in batch:
                                attributes_g=get_attributes(g)
                
                                if (attributes_f["read"]=="1"):

                                        if (attributes_f["ID"]==attributes_g["ID"]
                                        and attributes_f["status"]==attributes_g["status"]
                                        and attributes_g["read"]=="2"):
                                                pair=[f,g]
                                                paired.insert(0,pair)
                                                pair_count=pair_count+1
        print("Number of paired fastqs: "+str(pair_count))
    
        update_samplesheet(paired,all_csv)
               
                
        
        
        
        #update_samplesheet(batch,batch_name)
        
def get_attributes(samplesheet_name):
        attributes={}
        attributes["ID"]=re.search("P.+_L",samplesheet_name).group(0)[1:-3]
        attributes["sex"]=re.search("P.+_L",samplesheet_name).group(0)[-3:-2]

        if not (re.search("_clean",samplesheet_name)==None):

                attributes["read"]=re.search("_clean.+.fastq.gz",samplesheet_name).group(0)[6:7]
        else:
                attributes["read"]=re.search("_R._B",samplesheet_name).group(0)[2:3]
             
        attributes["batch"]=re.search("_B.+_tumis",samplesheet_name).group(0)[2:5]
        attributes["lane"]=re.search("_L.+_tumis.fastq",samplesheet_name).group(0)[4:5]
     

        if ((re.search("tumour",samplesheet_name)!=None)):
                attributes["status"]="tumour"
        elif ((re.search("normal",samplesheet_name)!=None)):
                attributes["status"]="normal"
 
        return attributes

def check_csv_for_text(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if any(cell.strip() for cell in line.split(',')):
                print(file_path)
                return True
    return False


def sort_samplesheet(samplesheet_name):
        
        data = pd.read_csv(samplesheet_name)
        sorted_data = data.sort_values('patient')
        sorted_data.to_csv(samplesheet_name, index=False)

def update_samplesheet(fastq_pairs, samplesheet_name):
        f=open(samplesheet_name,"a")
     
        writer=csv.writer(f, delimiter="|") 
        if not (check_csv_for_text(samplesheet_name)):
                writer.writerow(["patient,sex,status,sample,lane,fastq_1,fastq_2"])

        num_pairs=len(fastq_pairs)
      

        i=0
        for i in range(len(fastq_pairs)):
                ID=fastq_pairs[i][0]
 
                fastq_pair=fastq_pairs[i]
                if (("clean1.fastq.gz" in fastq_pair[0])):
                        read1=fastq_pair[0]
                        read2=fastq_pair[1]
                elif ("_R1_" in fastq_pair[0]):
                        read1=fastq_pair[0]
                        read2=fastq_pair[1]
                else:
                        read1=fastq_pair[1]
                        read2=fastq_pair[0]
                attributes_1=get_attributes(read1)
                attributes_2=get_attributes(read2)
                patient=attributes_1["ID"]
              
                if (attributes_1["lane"]==attributes_2["lane"]):
                        lane=attributes_1["lane"]
                else:
                        lane="ERROR"
                        print("ERROR with lane attribute")
                        return 0

                sex=attributes_1["sex"] 
                if (sex=="M"):
                        sex="XY"
                elif (sex=="F"):
                        sex="XX"
                if ("tumour" in attributes_1["status"]):
                        status="1"
                elif ("normal" in attributes_1["status"]):
                        status="0"
                   
                else:
                        print("AAAHHH")
                        return 0
                

                if ("normal2" in attributes_1):
                        sample=patient+"_"+attributes_1["status"]+"sample_2"
                elif("tumour2" in attributes_1):
                        sample=patient+"_"+attributes_1["status"]+"sample_2"
                else:   
                        sample=patient+"_"+attributes_1["status"]+"_sample" 
                
                if ("clean" in read1):
                        basedir=bbsplit_dir
                else:
                        basedir="/data/local/MD_scholarly/data/processed/exome/02_renamed"
                row=str(attributes_1["ID"]+","+sex+","+status+","+sample+","+lane+","+basedir+"/"+read1+","+basedir+"/"+read2)
                i=i+1
                writer.writerow([row])
         
                #writer.writerow([str(element+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[0]+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[1]+",auto")])
   

main()