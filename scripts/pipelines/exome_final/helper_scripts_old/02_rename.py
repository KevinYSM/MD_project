import pandas as pd
import os
import re
import shutil

trimmed_umis_folder="/home/ubuntu/data/local/MD_scholarly/data/processed/exome/01_trimmed_umis"
rename_table="/home/ubuntu/data/local/MD_scholarly/data/raw/exome/all/file_rename_tumis_table.csv"
sex_table="/data/local/MD_scholarly/data/raw/exome/all/patient_sex.csv"
#failed_list=["P126_L002_tumour1_R1_B952_tumis.fastq.gz","P126_L002_tumour1_R2_B952_tumis.fastq.gz","P37_L002_tumour1_R1_B693_tumis.fastq.gz","P37_L002_tumour1_R2_B693_tumis.fastq.gz"
  #           ,"P37_L002_normal1_R1_B693_tumis.fastq.gz","P37_L002_normal1_R2_B693_tumis.fastq.gz","P39_L001_tumour1_R1_B374_tumis.fastq.gz","P39_L001_tumour1_R2_B374_tumis.fastq.gz","P50_L001_tumour1_R1_B374_tumis.fastq.gz"
 #            ,"P50_L001_tumour1_R2_B374_tumis.fastq.gz","P54_L002_normal1_R1_B952_tumis.fastq.gz","P54_L002_normal1_R2_B952_tumis.fastq.gz"]

def verify_rename(old_file_loc, new_file_loc):
    
        #Verify Patient ID
        new_file_id=re.search("P.+_L",new_file_loc).group(0)[1:-3]
        if (re.search(new_file_id,old_file_loc)==None):
                print("Fail. Incorrect Patient ID: "+new_file_id)
                return False
        #Verify Lane
        old_file_lane=re.search("L0.+_R",old_file_loc).group(0)[:-2]
        if (re.search(old_file_lane,new_file_loc)==None):
                print("Fail. Incorrect lane: "+ old_file_lane)
                return False
        
        #Verify Sample. PDX, PBMC, BIOPSY, TIL
        if ((re.search("PDX",old_file_loc)!=None) or (re.search("biopsy",old_file_loc)!=None)):
                old_sample_type="tumour"
        elif ((re.search("TIL",old_file_loc)!=None) or (re.search("PBMC", old_file_loc)!=None)):
                old_sample_type="normal"
        else:
                print("Fail. Could not identify tumour/normal: "+old_file_loc)
                return False
        if (re.search(old_sample_type,new_file_loc)==None):
                print("######")
                print(old_file_loc)
                print(new_file_loc)
                print("Fail. Inorrect sample type: "+old_sample_type)
                
                return False
    
        #Verify Read
        old_file_read=re.search("_R._",old_file_loc).group(0)[1:3]
        if (re.search(old_file_read,new_file_loc)==None):
                print("Fail. Inorrect read: "+old_file_read)
                return False

        #Verify Batch
        #old_file_folder=re.search("GWA-JN-...",old_file_loc).group(0)[-3:]
        #if (re.search(old_file_folder,new_file_loc)==None):
        #        print("Fail. Could not identify batch: "+ old_file_folder)
        #        return False

        return True

def find_tumi(raw_name, tumi_dir):
        regex_search_string=raw_name[:-9]
        trimmed_umi_files=os.listdir(tumi_dir)

        for file in trimmed_umi_files:
                if not (re.search(regex_search_string,file)==None):
                        return file
                
        print("Oops:"+raw_name)

def main(trimmed_umis_folder):
        df=pd.read_csv(rename_table,names=["Patient_ID","Lane","Sample","Read","Batch","Raw_name","New_name","Tumi_name"])
        sex_df=pd.read_csv(sex_table, names=["Patient_ID","Sex"])
        sex_df["Patient_ID"]=pd.to_numeric(sex_df["Patient_ID"].str.strip())
        sex_df["Sex"]=sex_df["Sex"].str.strip()
        count=0
        print(len(df))
        for row in df.iterrows():
                
                raw_name=row[1]["Raw_name"]
                new_name_without_sex=row[1]["Tumi_name"]
                patient_id=int(re.search("P.+_L",new_name_without_sex).group(0)[1:-2])
               
                filtered_sex = sex_df.loc[sex_df["Patient_ID"] == patient_id, "Sex"]
                if not filtered_sex.empty:
                        patient_sex = filtered_sex.iloc[0]
                else:
                        patient_sex="U"
                
                sex_index=re.search("P.+_L",new_name_without_sex).end()-2
        
                new_name_with_sex=new_name_without_sex[:sex_index]+str(patient_sex)+new_name_without_sex[sex_index:]
               
            
                batch=row[1]["Batch"][1:]
                

                #Find old file
                old_file_loc=trimmed_umis_folder+"/"+find_tumi(raw_name, trimmed_umis_folder)
             

                #Make sure file exists 
                if not (os.path.isfile(old_file_loc)):
                        print("Oops!")
                        print(old_file_loc+" not found")
                        return
                else:
                        new_file_loc="/data/local/MD_scholarly/data/processed/exome/02_renamed/"+new_name_with_sex
                        
                        if (verify_rename(old_file_loc, new_file_loc)):
                              
                                count=count+1
                                if not (os.path.isfile(new_file_loc)):
                                        print("Was missing: "+new_file_loc)
                                        shutil.copyfile(old_file_loc, new_file_loc)
                                #print(new_file_loc)
                                        
        print(count)
        return 0
main(trimmed_umis_folder)

