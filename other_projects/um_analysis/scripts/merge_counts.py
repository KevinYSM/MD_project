import os
import pandas as pd
import re

raw_counts_folder="/data/local/MD_project/other_projects/um_analysis/data/raw"
def main():
        counts_files=os.listdir(raw_counts_folder)
        merged_df=read_file(os.path.join(raw_counts_folder, counts_files[0]))
        counts_files.pop(0)

        for file in counts_files:
                path=os.path.join(raw_counts_folder, file)
                new_df=read_file(path)
                merged_df.merge(new_df, )
        
        return 0

def read_file(file):
        sample_id=
        colnames=["gene_id",sample_id]
        return pd.read_csv(file,names=colnames)



def test():
        counts_files=os.listdir(raw_counts_folder)
        path=os.path.join(raw_counts_folder, counts_files[0])
        pandas_df=read_file(path)

        print(pandas_df)
        

main()