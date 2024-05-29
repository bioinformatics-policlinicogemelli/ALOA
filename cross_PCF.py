import pandas as pd 
import itertools
import json
from loguru import logger
import warnings
import os
import glob
import pathlib
from image_proc_functions import pheno_filt
import pandas as pd

pd.set_option('mode.chained_assignment', None)
warnings.simplefilter(action='ignore', category=FutureWarning)

def add_celltype(df, celltype_list):
    pheno_df=pheno_filt(df)
    pheno_df=pheno_df[['Cell.ID', 'Cell.X.Position', 'Cell.Y.Position', 'Pheno']]
    pheno_df=pheno_df.loc[pheno_df['Pheno'].isin(celltype_list["Phenotype"])]
    
    cell_type=celltype_list.set_index('Phenotype').to_dict('dict')['Cell_Type']
    pheno_df['Celltype']=pheno_df['Pheno'].replace(cell_type)
    
    return pheno_df
    
    
    
def main():
    print("\n################################## CROSS PCF ##################################\n")
    
    logger.info("Start cross PCF analysis: This step will provide the cross PCF evaluation for each ROI\n")
    
    with open("config.json") as f:
        data=json.load(f)
    
    input_folder = os.path.join(data["Paths"]["data_input_folder"],"raw_data")
    
    pathlib.Path(data["Paths"]["output_folder"]).mkdir(parents=True, exist_ok=True)
    output_folder = os.path.join(data["Paths"]["output_folder"],"Cross_PCF")
    pathlib.Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    celltype_list=pd.read_csv(os.path.join(data["Paths"]["data_input_folder"], "cellType_dict.tsv"), sep="\t")
    data_list=pd.read_csv(os.path.join(data["Paths"]["data_input_folder"], "sample_sheet.tsv"), sep="\t")
    groups=data_list["Group"].unique()
    
    for g in groups:
        logger.info(f"{len(groups)} found!")
        g_output=os.path.join(output_folder, g)
        pathlib.Path(g_output).mkdir(parents=True, exist_ok=True)
        
        data_g=data_list[data_list["Group"].isin([g])]
        
        for pzt in data_g["sbj_ID"]:
            pzt_path=os.path.join(input_folder, pzt)
            pzt_roi=[f for f in os.listdir(pzt_path) if f.endswith('cell_seg_data.txt')]
            pzt_roi=list(map(lambda x: os.path.join(pzt_path,x), pzt_roi))
            
            for ff in pzt_roi:
                df=pd.read_csv(ff, sep="\t")

                #reorganize DB
                pheno_df=add_celltype(df, celltype_list)
                
                combinations=list(itertools.combinations(celltype_list["Cell_Type"], 2))
                
                print("")

if __name__ == '__main__':
    main()