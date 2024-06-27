import os
import pandas as pd
import re

raw_counts_folder="/data/local/MD_project/other_projects/um_analysis/data/raw_renamed"
def main():
        counts_files=os.listdir(raw_counts_folder)
        merged_df=read_file(counts_files[0])
        counts_files.pop(0)

        for file in counts_files:
                
                new_df=read_file(file)
                merged_df=pd.merge(new_df, merged_df,on="gene_id", how="outer")
        merged_df.to_csv("merged_df.csv", index=False)
        return 0

def read_file(file_name):
        sample_id=file_name.split("-")[0]
        print(sample_id)
        path=os.path.join(raw_counts_folder, file_name)
        colnames=["gene_id",sample_id]

        df=pd.read_csv(path,names=colnames, sep="\t")
        
        return df


def test():
        counts_files=os.listdir(raw_counts_folder)
        path=os.path.join(raw_counts_folder, counts_files[0])
        pandas_df=read_file(path)

        print(pandas_df)
        

main()