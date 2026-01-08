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

from loguru import logger
import pathlib
import os
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math
import tap
import glob

#******************************************************
def raw_count_cells(PATH_MERGE_FOLDER,list_pheno):
    logger.info("Start raw count calculation:")
    dict_global_info = {}
    for _d in [f for f in os.listdir(PATH_MERGE_FOLDER) if not f.startswith('.')]:
        sub_directory=os.path.join(PATH_MERGE_FOLDER,_d)
        if not os.listdir(sub_directory):
            logger.error(f"Directory {_d} is empty")
            exit()
        for f in [f for f in os.listdir(sub_directory) if not f.startswith('.')]:
            _id_paz=f.replace("Merge_cell_seg_data_clean_","").replace(".txt","")
            logger.info(f"Subject: {_id_paz}")
            dict_global_info.setdefault(_d,{})
            dict_global_info[_d].setdefault(_id_paz,{})
            _file=os.path.join(sub_directory,f)
            _data=pd.read_csv(_file, sep="\t")
            _data["Pheno"]=list(map(lambda x: x.replace("+","+,"),_data["Pheno"]))
            _data["Pheno"]=list(map(lambda x: x[:-1],_data["Pheno"]))
            logger.info(f"Reading File {_file}")
            _total_cells=len(_data)
            dict_global_info[_d][_id_paz]["Total_Cells"]=_total_cells
            _data_groupped = _data.groupby(["Pheno"])["Pheno"].count()
            _data_dict = {frozenset(k.replace("OTHER,", "").replace(",OTHER", "").split(",")): v for k, v in _data_groupped.to_dict().items()}
            for main_pheno in list_pheno:
                dict_global_info[_d][_id_paz][f"{main_pheno}"] = _data_dict.get(frozenset(main_pheno.split(",")), 0)
    logger.info("End raw count calculation")
    return dict_global_info 

#******************************************************
def calculate_mean_group_cells(dictionary_raw_count):
    count_mean_group_cells={}
    for _k,_v in dictionary_raw_count.items():
        _temp=[]
        count_mean_group_cells.setdefault(_k,0)
        for _,_v2 in _v.items():
            _temp.append(_v2["Total_Cells"])
        count_mean_group_cells[_k]=round(np.mean(_temp),2)
    return count_mean_group_cells

#******************************************************
def normalized_count_cells(dictionary_raw_count,dictionary_mean_count):
    logger.info("Start normalized count calculation")
    dict_normalized_count_cells={}
    for group, patients in dictionary_raw_count.items():
        dict_normalized_count_cells[group]={}
        mean_group=dictionary_mean_count[group]
        for patient, counters in patients.items():
            logger.info(f"Subject: {patient}")
            dict_normalized_count_cells[group][patient]={}
            total_cells=counters["Total_Cells"]
            dict_normalized_count_cells[group][patient]["Total_Cells"]=total_cells
            for pheno, value in counters.items():
                if pheno=="Total_Cells": continue
                _new_value=(value/total_cells)*mean_group
                dict_normalized_count_cells[group][patient][pheno]=_new_value
    logger.info("End normalized count calculation")
    return dict_normalized_count_cells

#******************************************************
def create_output_dir(path_output,groups_name):
    for _g in groups_name:
        temp_folder=os.path.join(path_output,_g)
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
            logger.info(f" Created folder {temp_folder}")

#******************************************************
def create_summury_file(path_output_results,dictionary_count,type_data):
    logger.info("writing csv file with count info for each patient")
    for group,patients in dictionary_count.items():
        logger.info(f"Group: {group}")
        dire = os.path.join(path_output_results,group,"csv")
        pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        for patient,phenos in patients.items():
            logger.info(f"Subject: {patient}")
            with open(os.path.join(dire, f"{type_data}_count_{patient}.csv"),"w") as f:
                f.write(f"Patient\tPheno\tCount_{type_data}\n")
                logger.info(f"Start Writing csv file...")
                for p in phenos:
                    f.write(f'{patient}\t{p}\t{phenos[p]}\n')
                logger.info(f"End Writing {type_data}_count_{patient}.csv for {patient}")

#******************************************************
def create_norm_all_file(path_output_results,dictionary_norm_count):
    logger.info("Start normalized count calculatin on all groups")
    for group,patients in dictionary_norm_count.items():
        dire = os.path.join(path_output_results,group,"csv")
        pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        for patient,phenos in patients.items():
            with open(os.path.join(dire, "all_norm_count_"+patient+".csv"),'w') as f:
                f.write(f"Patient\tPheno\tCount_all_Norm\n")
                logger.info(f"Writing normalized count file on all groups for patient {patient}")
                for p in phenos:
                    f.write(f'{patient}\t{p}\t{phenos[p]}\n')
    logger.info("End normalized count calculation on all groups")

#******************************************************
def merge_counts(path_output_results, pattern, output_filename, phenos_list):
    all_files = []
    for group_dir in [d for d in os.listdir(path_output_results) if not d.startswith('.')]:
        csv_dir = os.path.join(path_output_results, group_dir, "csv")
        files = glob.glob(os.path.join(csv_dir, pattern))
        all_files.extend(files)

    if not all_files:
        logger.warning(f"No files found for pattern {pattern}")
        return

    # Load all files
# Load all files and reconstruct Group from folder structure
    df_list = []
    for f in all_files:
        group = pathlib.Path(f).parents[2].name
        tmp = pd.read_csv(f, sep="\t")
        tmp["Group"] = group
        df_list.append(tmp)

    df = pd.concat(df_list, ignore_index=True)


    # Infer the count column automatically
    count_col = [c for c in df.columns if c.startswith("Count_") or c.startswith("Count_all")][0]

    # Pivot to WIDE format
    df_wide = df.pivot_table(
        index=["Group", "Patient"],
        columns="Pheno",
        values=count_col,
        aggfunc="first"
    ).reset_index()

    # Ensure all phenotypes exist as columns (even if missing)
    for ph in phenos_list:
        if ph not in df_wide.columns:
            df_wide[ph] = 0

    if "Total_Cells" not in df_wide.columns:
        df_wide["Total_Cells"] = 0

    # Reorder columns
    ordered_cols = ["Group", "Patient"] + phenos_list + ["Total_Cells"]
    existing = [c for c in ordered_cols if c in df_wide.columns]
    df_wide = df_wide[existing]

    # Save final wide file
    merged_path = os.path.join(path_output_results, output_filename)
    df_wide.to_csv(merged_path, sep="\t", index=False)

    logger.info(f"✅ WIDE merged file saved as {merged_path}")

#******************************************************
def create_and_merge_counts(path_output_results, dict_raw_count, phenos_list,
                            dict_norm_count=None, dict_norm_all_count=None):

    create_summury_file(path_output_results, dict_raw_count, "Raw")
    merge_counts(path_output_results, "Raw_count_*.csv", "merged_raw.csv", phenos_list)

    if dict_norm_count is not None:
        create_summury_file(path_output_results, dict_norm_count, "Norm")
        merge_counts(path_output_results, "Norm_count_*.csv", "merged_norm.csv", phenos_list)

    if dict_norm_all_count is not None:
        create_norm_all_file(path_output_results, dict_norm_all_count)
        merge_counts(path_output_results, "all_norm_count_*.csv", "merged_all_norm.csv", phenos_list)

#******************************************************
def bar_plot(path_output_result,dict_data,type_data):
    logger.info(f"Barplot creation for {type_data} count")
    for group,data in dict_data.items():
        dire = os.path.join(path_output_result,group,"Bar_plot")
        pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        logger.info(f"Group: {group}")
        max_val_y=0
        fig=go.Figure()
        for patient,targets in data.items():
            logger.info(f"Subject: {patient}")
            sorted_target=dict(sorted(targets.items(),key=lambda x:x[0]))
            target=list(sorted_target.keys())
            target.remove("Total_Cells")
            values=list(sorted_target.values())[:-1]
            max_val=max(values)
            if max_val>max_val_y:
                max_val_y=max_val
            fig.add_trace(go.Bar(x=target,y=values, name=patient))
        fig.update_layout(barmode='group',title_text=f"Phenotypes {type_data}-{group}",xaxis=dict(
            title="Phenotypes",showticklabels=True),yaxis=dict(title=f'log({type_data} Count)',showticklabels=True,type="log"),legend_title_text="Patients")
        fig.update_yaxes(range=[0,(math.log(max_val_y,10)+1.0)])
        fig.write_image(os.path.join(dire,"Bar_Plot_"+type_data+".jpeg"),scale=6)

#******************************************************
def calculate_mean_total_groups(dictionary_raw_count):
    _count_total=[]
    for _k,_v in dictionary_raw_count.items():
        for _k2,_v2 in _v.items():
            _count_total.append(_v2["Total_Cells"])
    mean_cells=np.mean(_count_total)
    return mean_cells

#******************************************************
def normalized_count_on_all_groups(dictionary_raw_count,mean_all_groups):
    logger.info("Starting normalized count on all groups")
    dict_normalized_count_cells={}
    for group, patients in dictionary_raw_count.items():
        logger.info(f"Group: {group}")
        dict_normalized_count_cells[group]={}
        for patient, counters in patients.items():
            logger.info(f"Subject: {patient}")
            dict_normalized_count_cells[group][patient]={}
            total_cells=counters["Total_Cells"]
            dict_normalized_count_cells[group][patient]["Total_Cells"]=total_cells
            for pheno, value in counters.items():
                if pheno == "Total_Cells": continue
                _new_value=(value/total_cells)*mean_all_groups
                dict_normalized_count_cells[group][patient][pheno]=_new_value
    logger.info("End normalized count on all groups")
    return dict_normalized_count_cells

#******************************************************
def prepare_data_box_plot(dictionary_count):
    _data_all=pd.DataFrame()
    _data_norm=[]
    for _k,_v in dictionary_count.items():
        for _,_v2 in _v.items():
            for _t in _v2.keys():
                if _t=="Total_Cells": continue
                _data_norm.append([_k,_t,_v2[_t]])
    _data_all=pd.DataFrame(_data_norm,columns=["group","pheno","Count"])
    return _data_all

#******************************************************
def create_output_box_plot_dir(path_output):
    temp_folder=os.path.join(path_output,"Box Plots")
    pathlib.Path(temp_folder).mkdir(parents=True, exist_ok=True)

#******************************************************
def create_comparison_box_plot(path_output_result, data_all, p_adjust, test, type_data):
    logger.info(f"Box Plot creation on {type_data} count")
    x="pheno"
    y="Count"
    hue="group"
    hue_order=list(data_all["group"].unique())
    labels=sorted(list(data_all["pheno"].unique()))
    base_dir=os.path.join(path_output_result,"Box Plots",type_data)
    pathlib.Path(base_dir).mkdir(parents=True, exist_ok=True)
    # Mega plot
    filename=os.path.join(base_dir,f"box_plot_comparison_{type_data}.jpeg")
    _run_stat_test(data_all, labels, hue_order, p_adjust, test, x, y, hue, filename, type_data)
    # Single plots
    for pheno in labels:
        logger.info(f"Creating single boxplot for phenotype: {pheno}")
        df_single=data_all[data_all["pheno"]==pheno].copy()
        filename_single=os.path.join(base_dir,f"box_plot_{pheno}_{type_data}.jpeg")
        _run_stat_test(df_single, [pheno], hue_order, p_adjust, test, x, y, hue, filename_single, type_data)

#******************************************************
def _run_stat_test(data, labels, hue_order, p_adjust, test, x, y, hue, filename, type_data):
    kwargs_common={"width":4000,"height":1000,"title":f"Comparison {type_data} Count","log_y":True,
                   "labels":{"pheno":"Phenotypes","value":f"log({type_data} Counts)","group":"Group"}}
    try:
        if len(hue_order)==2 and test!="paired":
            logger.info(f"Mann-Whitney Test → {filename}")
            tap.plot_stats(data, x, y, order=labels, filename=filename,
                           type_correction=p_adjust, export_size=(1400,950,3),
                           subcategory=hue, kwargs=kwargs_common)
        elif len(hue_order)==2 and test=="paired":
            logger.info(f"Wilcoxon Test → {filename}")
            tap.plot_stats(data, x, y, order=labels, filename=filename,
                           type_test="wilcoxon", type_correction=p_adjust, export_size=(1400,950,3),
                           subcategory=hue, kwargs=kwargs_common)
        elif len(hue_order)>2:
            logger.info(f"Kruskal-Wallis Test → {filename}")
            if p_adjust is None: p_adjust="bonferroni"
            tap.plot_stats(data, x, y, type_test="dunn", order=labels,
                           subcategory=hue, filename=filename,
                           type_correction=p_adjust, export_size=(1400,950,3),
                           kwargs=kwargs_common)
        logger.info(f"✅ Saved: {filename}")
    except Exception as e:
        logger.error(f"❌ Error in statistical test for {filename}: {e}")

#******************************************************
def main(data):
    print("\n############################# DESCRIPTIVE ANALYSIS #############################\n")
    logger.info("Start overview process...")
    path_merge_folder=os.path.join(data["Paths"]["output_folder"],"Merged_clean")
    logger.info(f"Merge data found in {path_merge_folder}")
    if not os.listdir(path_merge_folder):
        logger.critical(f"{path_merge_folder} is an empty directory")
    path_output_results=os.path.join(data["Paths"]["output_folder"],"Descriptive")
    logger.info(f"Descriptive output will be stored in {path_output_results}")
    dict_pheno_interested=data["Phenotypes"]["pheno_list"]
    logger.info(f"The following phenotypes will be considered: {dict_pheno_interested}")
    dict_raw_count=raw_count_cells(path_merge_folder,dict_pheno_interested)
    groups=list(dict_raw_count.keys())
    logger.info(f"{len(groups)} group(s) found: {groups}")
    create_output_dir(path_output_results,groups)
    p_adjust=data["Stats"]["p_adj"].lower()
    if p_adjust=="": p_adjust=None
    test=data["Stats"]["sample_type"]

    # CREATE CSV and merge
    phenos_list = list(dict_pheno_interested)
    create_and_merge_counts(path_output_results, dict_raw_count, phenos_list)


    # Raw plots
    if data["Descriptive"]["raw"]:
        bar_plot(path_output_results, dict_raw_count, "Raw")
        if len(groups)>=2:
            create_output_box_plot_dir(path_output_results)
            df_raw=prepare_data_box_plot(dict_raw_count)
            create_comparison_box_plot(path_output_results, df_raw, p_adjust, test, "Raw")
        else:
            logger.critical("Found only one group - impossible to continue with groups comparison and box plot figure")

    # Normalized plots
    if data["Descriptive"]["normalized"]:
        mean_group=calculate_mean_group_cells(dict_raw_count)
        norm_count=normalized_count_cells(dict_raw_count, mean_group)
        create_summury_file(path_output_results, norm_count, "Norm")
        merge_counts(path_output_results, "Norm_count_*.csv", "merged_norm.csv", phenos_list)
        bar_plot(path_output_results, norm_count, "Normalized")
        if len(groups)>=2:
            mean_group_all=calculate_mean_total_groups(dict_raw_count)
            norm_count_all_groups=normalized_count_on_all_groups(dict_raw_count, mean_group_all)
            create_norm_all_file(path_output_results, norm_count_all_groups)
            merge_counts(path_output_results, "all_norm_count_*.csv", "merged_all_norm.csv", phenos_list)
            create_output_box_plot_dir(path_output_results)
            df_norm=prepare_data_box_plot(norm_count_all_groups)
            create_comparison_box_plot(path_output_results, df_norm, p_adjust, test, "Normalized")
        else:
            logger.critical("Found only one group-impossible to continue with group comparison")

    logger.info("End descriptive analysis!\n")
    return()

if __name__=="__main__":
    main()
