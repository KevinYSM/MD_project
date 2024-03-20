import pandas as pd
import os
import re
import sys

raw_exome_folder=sys.argv[1]
rename_table="/home/ubuntu/scratch/MD_project/scripts/pipelines/exome/file_rename_table.csv"
sex_table="/home/ubuntu/scratch/MD_project/scripts/pipelines/exome/patient_sex.csv"

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
        old_file_folder=re.search("GWA-JN-...",old_file_loc).group(0)[-3:]
        if (re.search(old_file_folder,new_file_loc)==None):
                print("Fail. Could not identify batch: "+ old_file_folder)
                return False

        return True

def main(raw_exome_folder):
        df=pd.read_csv(rename_table,names=["Patient_ID","Lane","Sample","Read","Batch","Old_name","New_name"])
        sex_df=pd.read_csv(sex_table, names=["Patient_ID","Sex"])
        sex_df["Patient_ID"]=pd.to_numeric(sex_df["Patient_ID"].str.strip())
        sex_df["Sex"]=sex_df["Sex"].str.strip()
        for row in df.iterrows():
                old_name=row[1]["Old_name"]
                new_name_without_sex=row[1]["New_name"]
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
                old_file_loc=raw_exome_folder+"/GWA-JN-"+batch+"/fastqs/"+old_name

                #Make sure file exists
                if not (os.path.isfile(old_file_loc)):
                        print("Oops!")
                        print(old_file_loc+" not found")
                        return
                else:
                        new_file_loc="/home/ubuntu/scratch/MD_project/data/exome/raw/unbatched/"+new_name_with_sex
                        
                        if (verify_rename(old_file_loc, new_file_loc)):
                                if (os.path.isfile(new_file_loc)):
                                        os.remove(new_file_loc)
                                os.symlink(old_file_loc, new_file_loc)
        return 0
main(raw_exome_folder)

