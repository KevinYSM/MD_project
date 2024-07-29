import pandas as pd
import os
import re
import shutil
import csv
import sys

bbsplit_dir="/home/ubuntu/data/local/MD_project/data/exome/processed/02_disambiguated_new"
trimmed_umis_dir=sys.argv[1]
batch_ids=["335","374","693","952"]
rename_table="/home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/helper_scripts/file_rename_tumis_table.csv"
sex_table="/home/ubuntu/data/local/MD_project/scripts/pipelines/exome/patient_sex.csv"
def main():
        #normal_fastqs=get_normal_fastqs(trimmed_umis_dir)
       # print(normal_fastqs[1][1])
        


        pattern=re.compile(r'.*\.fastq\.gz')
        all_fastqs=get_all_fastqs(trimmed_umis_dir)
        print(all_fastqs)
        create_all_samplesheet(all_fastqs)
       
        
        return 0



def get_batches(filtered_files, batch_ids):
        batches=[]

        for batch in batch_ids:
                filter=re.compile(re.escape("B"+batch))
                batch_files=[file for file in filtered_files if filter.search(file)]
                batches.append(batch_files)
        return batches

def get_tumour_fastqs(dir):
        df=pd.read_csv(rename_table,names=["Patient_ID","Lane","Sample","Read","Batch","Raw_name","New_name","Tumi_name"])
        sex_df=pd.read_csv(sex_table, names=["Patient_ID","Sex"])
        

        
        count=0
        tumour_fastqs_with_key=[]
        for basename in os.listdir(dir):
                f=os.path.join(dir, basename)
                if os.path.isfile(f) and f.endswith(".fastq.gz"):
                        
                        row=df.loc[df['Raw_name'].str.contains(basename[:-29])]
                       
                        if ("tumour" in str(row["Sample"].values)):                                
                                tumour_fastqs_with_key.append([basename,row["Tumi_name"].values.tolist()[0]])
                                count=count+1
                                
                                
                                
       
        return tumour_fastqs_with_key

def get_all_fastqs(dir):
        df=pd.read_csv(rename_table,names=["Patient_ID","Lane","Sample","Read","Batch","Raw_name","New_name","Tumi_name"])
        sex_df=pd.read_csv(sex_table, names=["Patient_ID","Sex"])
        

        
        count=0
        failed_fastqs=[]
        all_fastqs=[]
        for basename in os.listdir(dir):
                f=os.path.join(dir, basename)
                if os.path.isfile(f) and f.endswith(".fastq.gz"):
                        print(basename)
                        print(basename[:-29])
                     
                        row=df.loc[df['Raw_name'].str.contains(basename[:-29])]     
                                
                       
                        

                        if not (row.empty):
                                tumi_name=row["Tumi_name"].values.tolist()[0]
                                all_fastqs.append(tumi_name)
                                count=count+1
                        else:
                                failed_fastqs.append(basename)
                        
                        
                                  
        with open("/home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/helper_scripts/failed_fastqs.csv", 'w') as myfile:
                wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
                for fastq in failed_fastqs:
                        wr.writerow([fastq])
        return all_fastqs

def get_normal_fastqs(dir):
        df=pd.read_csv(rename_table,names=["Patient_ID","Lane","Sample","Read","Batch","Raw_name","New_name","Tumi_name"])
        sex_df=pd.read_csv(sex_table, names=["Patient_ID","Sex"])
        

        
        count=0
        normal_fastqs_with_key=[]
        for basename in os.listdir(dir):
                f=os.path.join(dir, basename)
                if os.path.isfile(f) and f.endswith(".fastq.gz"):
                        
                        row=df.loc[df['Raw_name'].str.contains(basename[:-29])]
                       
                        if ("normal" in str(row["Sample"].values)):                                
                                normal_fastqs_with_key.append([basename,row["Tumi_name"].values.tolist()[0]])
                                count=count+1
                                
                                
                                
       
        return normal_fastqs_with_key


def get_fastqs(trimmed_umis_dir):
        
        tumour_samples_paired= get_tumour_fastqs(trimmed_umis_dir)
        tumour_samples_tumi_name=[]

        
        normal_samples_paired= get_normal_fastqs(trimmed_umis_dir)
        normal_samples_tumi_name=[]
        for f in normal_samples_paired:
                normal_samples_tumi_name.append(f[1])
        for f in tumour_samples_paired:
                tumour_samples_tumi_name.append(f[1])
        all_samples=tumour_samples_tumi_name+normal_samples_tumi_name

        
        print(len(tumour_samples_tumi_name),len(normal_samples_tumi_name))
        filter=re.compile(r'\.fastq\.gz$')
       
        return all_samples
        


       # all_csv="/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/all.csv"
        #if (os.path.exists(all_csv)):
        #        f=open(all_csv,"w")
       #         f.truncate()
        #        f.close()

        #batches=get_batches(filtered_files,batch_ids)
        
        #generate_samplesheets(batches)

        #sort_samplesheet("/data/local/MD_scholarly/scripts/pipelines/exome/04_samplesheets/samplesheets/all.csv")

        return 0

           


def create_all_samplesheet(all_samples):
        print(len(all_samples))
        all_csv="/home/ubuntu/data/local/MD_project/scripts/pipelines/exome_final/steps/02_disambiguate/tumi_tumour.csv"
        if (os.path.exists(all_csv)):
                f=open(all_csv,"w")
                f.truncate()
                f.close()
        paired=[]
        pair_count=0
        for sample_1 in all_samples:
                batch_id=get_attributes(sample_1)["batch"]
                
               
               
                attributes_1=get_attributes(sample_1)
                for sample_2 in all_samples:
                        attributes_2=get_attributes(sample_2)
        
                        if (attributes_1["read"]=="1"):
                            
                                if (attributes_1["ID"]==attributes_2["ID"]
                                and attributes_1["status"]==attributes_2["status"]
                                and attributes_2["read"]=="2"
                                and attributes_1["lane"]==attributes_2["lane"]):
                                        
                                        pair=[sample_1,sample_2]
                                        paired.insert(0,pair)
                                        pair_count=pair_count+1
                        
                                        
        print("Number of paired fastqs: "+str(pair_count))
    
        update_samplesheet(paired,all_csv)
               
                
        
        
        
        #update_samplesheet(batch,batch_name)
        
def get_sex(patient_id):

      
        sex_df=pd.read_csv(sex_table, names=["Patient_ID","Sex"])
        sex_df["Patient_ID"]=sex_df["Patient_ID"].str.strip()
        sex_df["Sex"]=sex_df["Sex"].str.strip()
 
        return sex_df.loc[sex_df["Patient_ID"]==patient_id]["Sex"].tolist()[0]

  
def get_attributes(samplesheet_name):
      
      
        attributes={}
     
        if (re.search("out_human_P.+_L",samplesheet_name)!=None):
                attributes["ID"]=re.search("out_human_P.+_L",samplesheet_name).group(0)[11:-2]
        else:
                attributes["ID"]=re.search("P.+_L",samplesheet_name).group(0)[1:-2]
     
        attributes["sex"]=get_sex(attributes["ID"])
        if  (re.search("out_human",samplesheet_name)!=None):
                attributes["read"]=re.search("_tumis_.+.fastq.gz",samplesheet_name).group(0)[7]
               
        else:
                attributes["read"]=re.search("_R._B",samplesheet_name).group(0)[2:3]
             
        attributes["batch"]=re.search("_B.+_tumis",samplesheet_name).group(0)[2:5]


        if (re.search("_L.+_tumis_1.fastq.gz",samplesheet_name)!=None):

                attributes["lane"]=re.search("_L.+_tumis_1.fastq.gz",samplesheet_name).group(0)[4:5]
        elif ((re.search("_L.+_tumis_2.fastq.gz",samplesheet_name)!=None)):
                attributes["lane"]=re.search("_L.+_tumis_2.fastq.gz",samplesheet_name).group(0)[4:5]
        else:
                attributes["lane"]=re.search("_L.+_tumis.fastq.gz",samplesheet_name).group(0)[4:5]
     

        if ((re.search("tumour",samplesheet_name)!=None)):
                attributes["status"]="tumour"
        elif ((re.search("normal",samplesheet_name)!=None)):
                attributes["status"]="normal"\
        
 
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
        normal_fastqs=get_normal_fastqs(trimmed_umis_dir)
        tumour_fastqs=get_tumour_fastqs(trimmed_umis_dir)
      
        writer=csv.writer(f, delimiter="|") 
        if not (check_csv_for_text(samplesheet_name)):
                writer.writerow(["patient,sex,status,sample,lane,fastq_1,fastq_2"])

        num_pairs=len(fastq_pairs)
      

        i=0
        for i in range(len(fastq_pairs)):
                ID=fastq_pairs[i][0]
                
                fastq_pair=fastq_pairs[i]
                
                if (("1.fastq.gz" in fastq_pair[0])):
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
                        print(read1,read2)
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
                
                if ("normal2" in read1):
                        sample=patient+"_"+attributes_1["status"]+"_sample_2"
                elif("tumour2" in read1):
                        sample=patient+"_"+attributes_1["status"]+"_sample_2"
                else:   
                        sample=patient+"_"+attributes_1["status"]+"_sample" 
                
                if ("normal" in read1):
                        for f in normal_fastqs:
                                tumi_name=f[1]
                       
                                if (tumi_name == read1):
                                        read1=trimmed_umis_dir+"/"+f[0]
                else:
                        for f in tumour_fastqs:
                                tumi_name=f[1]
                       
                                if (tumi_name==read1):
                                        read1=trimmed_umis_dir+"/"+f[0]
                                        print(read1)
                        
                if ("normal" in read2):
                        for f in normal_fastqs:
                                tumi_name=f[1]
                                if (tumi_name == read2):
                                        read2=trimmed_umis_dir+"/"+f[0]
                else:
                        for f in tumour_fastqs:
                                tumi_name=f[1]

                                
                                if (tumi_name == read2):
                                        read2=trimmed_umis_dir+"/"+f[0]
                                        print(read2)
                       
                        

                basedir="/data/local/MD_scholarly/data/processed/exome/02_renamed"
                
                row=str(attributes_1["ID"]+","+sex+","+status+","+sample+","+lane+","+read1+","+read2)
                
                i=i+1
                writer.writerow([row])
         
                #writer.writerow([str(element+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[0]+","+"/data/local/MD_scholarly/data/raw/rnaseq/fastqs/"+fastq_pair[1]+",auto")])
   

main()