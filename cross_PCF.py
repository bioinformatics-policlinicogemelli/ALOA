import pandas as pd 
import itertools
import json
from loguru import logger
import warnings
import os
import re
import pathlib
from pcf_functions import *

pd.set_option('mode.chained_assignment', None)
warnings.simplefilter(action='ignore', category=FutureWarning)
 
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
    
    logger.info(f"{len(groups)} Group(s) found!")

    for g in groups:
        logger.info(f"Analyzing group: {g}")
        g_output=os.path.join(output_folder, g)
        pathlib.Path(g_output).mkdir(parents=True, exist_ok=True)
        
        data_g=data_list[data_list["Group"].isin([g])]
        
        combinations=list(itertools.combinations(celltype_list["Cell_Type"], 2))
        logger.info(f"{len(combinations)} cell type combination(s) found: {combinations}")
        
        #set radius of interest
        radiusOfInterest = data["pcf"]["radiusOfInterest"]
        logger.info(f"Radius {radiusOfInterest} was selected. It is recommended to choose an r value at least equal to or greater than 2 or 3 times the diameter of the cells under consideration.")
        
        for C_1, C_2 in combinations:

            logger.info(f"Analyzing combination: {C_1} - {C_2}")
            output_csv= create_output_csv(os.path.join(g_output,'csv'), C_1.replace(" ","_"), C_2.replace(" ","_"))
            #comb_output=os.path.join(g_output,C_1.replace(" ","_")+"-"+C_2.replace(" ","_"))

            for pzt in data_g["sbj_ID"]:
                logger.info(f"Analyzing subject: {pzt}")

                pzt_output=os.path.join(g_output, pzt)
                pathlib.Path(pzt_output).mkdir(parents=True, exist_ok=True)
                logger.info(f"output subject folder created in {pzt_output} path")

                comb_output=os.path.join(pzt_output,C_1.replace(" ","_")+"-"+C_2.replace(" ","_"))
                pathlib.Path(comb_output).mkdir(parents=True, exist_ok=True)
                
                pzt_path=os.path.join(input_folder, pzt)
                pzt_roi=[f for f in os.listdir(pzt_path) if f.endswith('cell_seg_data.txt') and not f.startswith('Merge')]
                pzt_roi=list(map(lambda x: os.path.join(pzt_path,x), pzt_roi))
                
                for ff in pzt_roi:
                    
                    roi_name=re.search(r'\[(.*?)\]', ff).group(1)
                    logger.info(f"Analyzing ROI: {roi_name}")
                    df=pd.read_csv(ff, sep="\t")

                    #reorganize DB
                    pheno_df=add_celltype(df, celltype_list, data["pcf"]["pheno_list"])
                    
                    if C_1 not in pheno_df['Celltype'].values or C_2 not in pheno_df['Celltype'].values:
                        print(f"{C_1} o {C_2} cell types not found in 'Celltype' column of {ff}. Skip to next combination.")
                        continue

                    roi_output=os.path.join(comb_output, "ROI_"+roi_name)
                    pathlib.Path(roi_output).mkdir(parents=True, exist_ok=True)

                    pc = load_point_cloud(pheno_df, roi_output)
                    
                    if data["pcf"]["all_pcf"]:
                        pc = all_cross_pcf(pc, pzt_output, roi_name, data["pcf"]["maxR"], data["pcf"]["annulusStep"], data["pcf"]["annulusWidth"])
                    
                    pcf_value_at_radius =selected_PCF(C_1, C_2, pc, roi_output, radiusOfInterest, data["pcf"]["maxR"], data["pcf"]["annulusStep"], data["pcf"]["annulusWidth"])
                    
                    TCM(C_1, C_2, radiusOfInterest, pc, roi_output)

                    count_C1 = len(pheno_df[pheno_df['Celltype'] == C_1])
                    count_C2 = len(pheno_df[pheno_df['Celltype'] == C_2])

                    append_to_csv(output_csv, pzt, roi_name, pcf_value_at_radius, count_C1, count_C2)
    
    logger.info("End PCF evaluation!\n")

if __name__ == '__main__':
    main()