#!/usr/bin/env python3
import pandas as pd
import os
import re
import shutil
import sys

rename_table="/home/ubuntu/scratch/MD_project/scripts/pipelines/exome/file_rename_tumis_table.csv"


def return_batch(current_name):
    basename=os.path.basename(current_name)

    df=pd.read_csv(rename_table,names=["Patient_ID","Lane","Sample","Read","Batch","Raw_name","New_name","Tumi_name"])

    basename_substring=basename[0:re.search("_R._",basename).span()[1]].strip()

    batch=df[df["Raw_name"].str.contains(basename_substring.strip())]["Batch"].values[0]
    
  
    return batch[1:]
   


def main(args):
    #print(args[0])
    print(return_batch(args[0]))
    
    return (return_batch(args[0]))

if __name__ == "__main__":
    main(sys.argv[1:])