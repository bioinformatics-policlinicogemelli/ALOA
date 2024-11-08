#Copyright 2024 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

import pandas as pd 
import itertools
from loguru import logger
import warnings
import os
import re
import pathlib
from pcf_functions import *
from pandas.api.types import is_numeric_dtype

pd.set_option('mode.chained_assignment', None)
warnings.simplefilter(action='ignore', category=FutureWarning)

def main(data=[]):
    
    logger.info("\n################################## CROSS PCF ##################################\n")
    
    logger.info("Start cross PCF analysis: This step will provide the cross PCF evaluation for each ROI\n")
    
    input_folder = os.path.join(data["Paths"]["data_input_folder"],"raw_data")
    
    pathlib.Path(data["Paths"]["output_folder"]).mkdir(parents=True, exist_ok=True)
    output_folder = os.path.join(data["Paths"]["output_folder"],"Cross_PCF")
    pathlib.Path(output_folder).mkdir(parents=True, exist_ok=True)

    celltype_list=load_df(data["Paths"]["data_input_folder"], "cellType_dict.tsv")
    if celltype_list.isnull().values.any():
        logger.critical("One or more Nan found in df. Check cellType_dict.tsv. This error could be due to the presence of double tabs between occurence inside the tsv file")
        return()
    
    celltype_list["Phenotype"]=list(map(lambda x: x.replace(",",""),celltype_list["Phenotype"]))

    if not celltype_list.columns.isin(["Phenotype", "Cell_Type"]).all():
        logger.critical("It seems that one or more column names of cellType_dict.tsv is not correct. Check that the columns are named as 'Phenotype' and 'Cell_Type'")
        raise(Exception("Error while reading cellType_dict.tsv"))
    
    data_list=load_df(data["Paths"]["data_input_folder"], "sample_sheet.tsv")
    
    if not data_list.columns.isin(["sbj_ID", "Group"]).all():
        logger.critical("It seems that one or more column names of cellType_dict.tsv is not correct. Check that the columns are named as 'sbj_ID' and 'Group'")
        raise(Exception("Error while reading sample_sheet.tsv"))
    
    groups=data_list["Group"].unique()
    
    logger.info(f"{len(groups)} Group(s) found!")

    #set radius of interest
    radiusOfInterest = data["Cross_pcf"]["radiusOfInterest"]
    logger.info(f"Radius={radiusOfInterest} micron was selected. It is recommended to choose an r value at least equal to or greater than 2 or 3 times the diameter of the cells under consideration.")
    
    if not data["Cross_pcf"]["only_stat"]:
        
        for g in groups:
            logger.info(f"Analyzing group: {g}")
            g_output=os.path.join(output_folder, g)
            pathlib.Path(g_output).mkdir(parents=True, exist_ok=True)
            
            data_g=data_list[data_list["Group"].isin([g])]
            
            combinations=list(itertools.combinations(celltype_list["Cell_Type"], 2))
            
            logger.info(f"{len(combinations)} cell type combination(s) found: {combinations}")
            
            for C_1, C_2 in combinations:

                logger.info(f"Analyzing combination: {C_1} - {C_2}")
                output_csv= create_output_csv(os.path.join(g_output,'summary'), C_1.replace(" ","_"), C_2.replace(" ","_"))

                for pzt in data_g["sbj_ID"]:
                    logger.info(f"Analyzing subject: {pzt}")

                    pzt_output=os.path.join(g_output, pzt)
                    comb_output=os.path.join(pzt_output,C_1.replace(" ","_")+"-"+C_2.replace(" ","_"))
                    
                    if data["Cross_pcf"]["save_images"]:
                        pathlib.Path(pzt_output).mkdir(parents=True, exist_ok=True)
                        logger.info(f"output subject folder created in '{pzt_output}' path")
                        
                        pathlib.Path(comb_output).mkdir(parents=True, exist_ok=True)
                    
                    pzt_path=os.path.join(input_folder, pzt)
                    pzt_roi=[f for f in os.listdir(pzt_path) if f.endswith('cell_seg_data.txt') and not f.startswith('Merge')]

                    pzt_roi=list(map(lambda x: os.path.join(pzt_path,x), pzt_roi))
                    
                    for ff in pzt_roi:
                        
                        roi_name=re.search(r'\[(.*?)\]', ff).group(1)
                        #logger.info(f"Analyzing ROI: {roi_name}")

                        try:
                            df=pd.read_csv(ff, sep="\t")
                            df.columns=list(map(lambda x: x.replace(" ",".") ,df.columns))
                            if not is_numeric_dtype(df["Cell.X.Position"]):
                                df=pd.read_csv(ff, sep="\t", decimal=",")
                                df.columns=list(map(lambda x: x.replace(" ",".") ,df.columns))
                        except FileNotFoundError:
                            logger.warning(f"No cell_segmentation file found in {ff}. Skip to the next one!")
                            continue

                        #reorganize DB
                        pheno_df=add_celltype(df, celltype_list, data["Cross_pcf"]["pheno_list"])
                        
                        if C_1 not in pheno_df['Celltype'].values or C_2 not in pheno_df['Celltype'].values:
                            logger.warning(f"{C_1} or {C_2} cell types not found in 'Celltype' column of {ff}. Skip to next combination.")
                            continue
                        
                        roi_output=os.path.join(comb_output, "ROI_"+roi_name)
                        rad_folder=os.path.join(roi_output,"r_"+str(radiusOfInterest))
                        
                        if data["Cross_pcf"]["save_images"]:
                            pathlib.Path(roi_output).mkdir(parents=True, exist_ok=True)
                            pathlib.Path(rad_folder).mkdir(parents=True, exist_ok=True)
                        
                        pc = load_point_cloud(pheno_df, rad_folder, C_1, C_2, data["Cross_pcf"]["save_images"])
                        
                        if data["Cross_pcf"]["all_pcf"]:
                            try:
                                pc = all_cross_pcf(pc, pzt_output, roi_name, data["Cross_pcf"]["maxR"], data["Cross_pcf"]["annulusStep"], data["Cross_pcf"]["annulusWidth"], data["Cross_pcf"]["save_images"])
                            except RuntimeError:
                                logger.error("Negative areas calculated! Skip to the next ROI")
                                continue
                        try:
                            pcf_value_at_radius =selected_PCF(C_1, C_2, pc, rad_folder, radiusOfInterest, data["Cross_pcf"]["maxR"], data["Cross_pcf"]["annulusStep"], data["Cross_pcf"]["annulusWidth"],data["Cross_pcf"]["save_images"])
                        except RuntimeError:
                            logger.error("Negative areas calculated! Skip to the next ROI")
                            continue
                        
                        if data["Cross_pcf"]["save_images"]:
                            tcm=TCM(C_1, C_2, radiusOfInterest, pc, rad_folder)
                            
                            if data["Cross_pcf"]["on_roi"]:
                                input_roi = os.path.join(input_folder.replace("raw_data","img_match"), os.path.basename(ff).replace("cell_seg_data.txt","composite_image"))
                                
                                tcm_on_roi(input_roi, tcm, radiusOfInterest, rad_folder, C_1, C_2)

                        count_C1 = len(pheno_df[pheno_df['Celltype'] == C_1])
                        count_C2 = len(pheno_df[pheno_df['Celltype'] == C_2])

                        append_to_csv(output_csv, pzt, roi_name, pcf_value_at_radius, count_C1, count_C2)

    #statistical analysis
    combinations=list(itertools.combinations(celltype_list["Cell_Type"], 2))
    
    if len(groups)<=1:
        logger.info("No Statistical analysis will be conducted: Not enough groups!")
    
    else:
        stat_folder=os.path.join(output_folder,"stats")
        pathlib.Path(stat_folder).mkdir(parents=True, exist_ok=True)

        stats_file=os.path.join(stat_folder,"stat_analysis_r_"+str(radiusOfInterest)+".tsv")
        create_stats_file(groups, stats_file)

        for C_1, C_2 in combinations:

            logger.info ("\nStart statistical analysis for combination " + C_1 + " - " + C_2)
            df_comb=pd.DataFrame()
            df_counts=pd.DataFrame()
            results=[C_1, C_2]
            for g in groups:
                
                summary_path=os.path.join(output_folder, g, "summary")
                
                try:
                    sum_comb=pd.read_csv(os.path.join(summary_path, C_1.replace(" ","_")+"-"+C_2.replace(" ","_")+".tsv"), sep="\t")
                    df_comb[f"PCF_{g}"]=sum_comb["PCF_r"]
                    df_counts[f"C1_{g}"]=sum_comb["Counts_C1"]
                    df_counts[f"C2_{g}"]=sum_comb["Counts_C2"]
                except NameError:
                    logger.warning(f"No summary combination file found in group {g}!")
                    continue
                except FileNotFoundError:
                    name=os.path.join(summary_path, C_1.replace(" ","_")+"-"+C_2.replace(" ","_")+".tsv")
                    logger.error(f"No such file or directory {name}")
                    continue
                
            if df_comb.shape==(0,0):
                continue
            
            if len(df_comb.columns)>=2:
                try:
                    pval, pvals=stats_eval(df_comb, groups, data["Stats"]["sample_type"], data["Stats"]["p_adj"].lower())
                except ValueError:
                    logger.error(f"Something went wrong with statistical test for {C_1}-{C_2}")
                    pval, pvals=np.nan, ['nc']
                except TypeError:
                    logger.error(f"Something went wrong with statistical test for {C_1}-{C_2}: it can be that one or more columns are all Nan")
                    pval, pvals=np.nan, ['nc']
            
            fill_stats_file(results, pd.concat([df_counts, df_comb], axis=1), stats_file, groups, pval, pvals)
            

    logger.info("End PCF evaluation!\n")

if __name__ == '__main__':
    main()