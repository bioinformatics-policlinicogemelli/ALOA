from datetime import datetime
from loguru import logger
import sys
import pathlib
import os
import json
import numpy as np
import polars as pl
import math
import pandas as pd
import tap
from scipy import stats
import itertools
import plotly.figure_factory as ff

def standardization_distance_all_image(values,paz):
    '''
    function to calculate z-score from distance raw count (distance value-mean(distances)/standard deviation))
    Parameters
    ----
    values : array
    paz : str

    Return
    ----
    '''
    if len(values)<2:
        logger.info(f"only one or less distance for {paz}")
        #print(f"only one or less distance for {paz}")
        return [],0,0

    mean=np.mean(values)
    std=np.std(values)

    if std==0:
        logger.info(f"for {paz} std = 0")
        #print(f"for {paz} std = 0")
        return [],0,0
                        
    values_standard= list(map(lambda x: (x - mean)/std,values))
    return values_standard,mean,std



def prepare_dataframe_distances(root_folder,pheno_from,pheno_to):
    '''
    function to create a dictionary with zeta score values for each group
    PARAMETERS
    ----
    root_folder:str
    pheno_from:str
    pheno_to:str
    ----
    '''

    DICT_DISTANCES = dict()

    #list_groups = os.listdir(root_folder)
    for group in [f for f in os.listdir(root_folder) if not f.startswith('.')]:
        if group=="Stats":
            continue
        #print("Parsing Group:",group)
        DICT_DISTANCES[group] = []
        if not os.listdir(os.path.join(root_folder,group)):
            logger.error(f"{os.path.join(root_folder,group)} is an empty folder")
            exit()
        for file in [f for f in os.listdir(os.path.join(root_folder,group)) if not f.startswith('.')]:
        #for directory in os.listdir(os.path.join(root_folder,group)):
            # if directory==".DS_Store":
                # continue
            # if directory=="csv":
            # if not os.listdir(os.path.join(root_folder, group,directory)):
            #     logger.error(f"{os.path.join(root_folder, group,directory)} is an empty folder")
            #     exit()
            # list_files=os.listdir(os.path.join(root_folder, group,directory))
            # for file in list_files:
            paz=file.replace("_Distance.txt","")
            file_path = os.path.join(root_folder, group, file)
            # file_path = os.path.join(root_folder, group,directory, file)

        #CALCOLA HEADER
            df = pl.scan_csv(file_path,separator="\t")
            df_columns = df.columns

            real_pheno_to = None
            set_pheno_to = set(pheno_to.split(","))
            for col in df_columns:
                set_col = set(col.replace("Distance to ", "").split(","))
                if set_col == set_pheno_to:
                    real_pheno_to = col
                    break

            if real_pheno_to is None:
                logger.warning(f"{pheno_to} not present as Distance_to for {paz}")
                continue

            #CERCO VALORE PHENO FROM
            real_pheno_from = None
            set_pheno_from=set(pheno_from.split(","))

            unique_values = pl.scan_csv(
                    file_path,
                    separator="\t",
                    null_values=["NA"]
                ).select(["Phenotype"]) \
                .drop_nulls() \
                .unique() \
                .collect()["Phenotype"] \
                .to_list()

            for val in unique_values:
                set_val = set(val.split(","))
                if set_val == set_pheno_from:
                    real_pheno_from = val
                    break

            if real_pheno_from is None:
                logger.warning(f"{pheno_from} not present in 'Phenotype' column for {paz}")
                continue

            #PRENDI DATI
            df_filtered = pl.scan_csv(
                    file_path,
                    separator="\t",
                    null_values=["NA"]
                ).select(["Phenotype", real_pheno_to]) \
                .drop_nulls() \
                .filter(pl.col("Phenotype") == real_pheno_from) \
                .select([real_pheno_to])
            
            
            values = list(map(lambda x: float(x),df_filtered.collect()[real_pheno_to].to_list()))
            values = [item for item in values if not(math.isnan(item)) == True]

            try:
                values_standard,mean,std= standardization_distance_all_image(values,paz)
                #if values_standard !=[]:
                    #f.write(f'{paz}\t{pheno_from}\t{pheno_to}\t{mean}\t{std}\n')
                DICT_DISTANCES[group] += values_standard
            except:
                continue

    return DICT_DISTANCES
#*****************************************************************
def create_output_dir(path_output_results):
    '''
    Funzione per la creazione dell cartella dove veranno inseriti gli output della descrittiva---> CARTELLA OUTPUT_DESCRIPTIVE con all'interno sottocartelle per ogni gruppo

    Parameters
    ----
    path_output : str
    groups_name : list

    Return
    ----
    None

    '''
 
    temp_folder=path_output_results
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

        #---->logger
        logger.info(f"Created folder {temp_folder}")
        
    else:
        #--->logger
       logger.info(f"Folder {temp_folder} already exists-deleting its contents for the new analysis")
       for direct in os.listdir(temp_folder):
            direct_path = os.path.join(temp_folder, direct)
            if not os.path.isfile(direct_path):
                for file in os.listdir(direct_path):
                    file_path=os.path.join(direct_path,file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                        #logger.info(f"File '{file_path}' deleted.")
                os.rmdir(direct_path)

#*****************************************************************
def create_df_distances(DICT_DISTANCES,path_output,pheno_from,pheno_to,save_zetascore=False):
    ##FIXME Metterlo come opzione se si vuole salvare il file csv
   

    pre_df = []
    for k, v in DICT_DISTANCES.items():
        for e in v:
            pre_df.append((k, e))

    df = pd.DataFrame(pre_df, columns=["GROUP", "DISTANCE"])
    if not os.path.exists(path_output):
        os.mkdir(path_output)
        df.to_csv(path_output, sep="\t",index=False)
    if len (df)!=0 and save_zetascore:
        direc=os.path.join(path_output,"csv")
        if not os.path.exists(direc):
            os.makedirs(direc)
            df.to_csv(f'{direc}/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv', sep="\t",index=False)
        df.to_csv(f'{direc}/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv', sep="\t",index=False)

    else:
        logger.warning(f"Dataframe is empty for distance from {pheno_from} to {pheno_to}")

    return df


def calculate_median_distribution(dictionary_group,groups):
    ##FIXME non salvare il file da solo con i valori della mediana
  
    dict_median={}
    
    for g in groups:
        #dict_meadian[g]=None
        dict_median[g]="NaN"
    
    for _k,_v in dictionary_group.items():
        if len (_v)==0:
            continue
        mean=np.median(_v)
        dict_median[_k]=mean
        
    grade_major=""

    #if any(value is None for value in dict_meadian.values()):
    if all(value=="NaN"for value in dict_median.values()):
        #print("both group's values are NaN for meadian distance calculation")
        grade_major="NaN"
        return grade_major,dict_median
    else:
        valid_value={k:v for k,v in dict_median.items() if str(v) !="NaN"}
        grade_major = max(valid_value, key=valid_value.get)
        return grade_major,dict_median


#*****************************************************************

def box_plots_distances(path_ouput_results,df,pheno_from,pheno_to):
   
    #dire=os.path.join(path_ouput_results,"box_plot")
    #if not os.path.exists(dire):
        #os.makedirs(dire)
    if len(df["GROUP"].unique())!=1: 
        dire=os.path.join(path_ouput_results,"box_plot")
        if not os.path.exists(dire):
            os.makedirs(dire)
        tap.plot_stats(df,x="GROUP",y="DISTANCE",filename=f'{dire}/distances_box_plot_{pheno_from}_to_{pheno_to}.png',kwargs={"title":f'Distance from {pheno_from} to {pheno_to}',"labels":{"GROUP": "Group",
                     "DISTANCE": r'$Distance_{z}$',
                     "GROUP": 'Group'
                 }})

#*****************************************************************
def statistical_test(path_output_result,df,pheno_from,pheno_to):
  
    groups=list(df["GROUP"].unique())

    #array contenente i valori delle distanze per le diverse cellulle per i diversi pazienti nei diversi gradi [[valori grado II],[valori grado III]]
    values_distance=[]

    for g in groups:
        temp=list(df[df["GROUP"]==g]["DISTANCE"].values)
        values_distance.append(temp)

    p_value=10

    if len(df["GROUP"].unique())==1:
        logger.warning("Only One Group - not statistical is possible")
    
    #Caso in cui abbiamo 2 GRUPPI---> Test di mannwhitney
    elif len(df["GROUP"].unique())==2:
        logger.info("Executive Mann-Whitney test")
        t_stat, p_value = stats.mannwhitneyu(values_distance[0], values_distance[1])
        #dict_p_value[f'{pheno_from}to{pheno_to}']=p_value
    
    #Caso in cui i gruppi sono più di due---> Test di Kruskal per più sample
    elif len(df["GROUP"].unique())>2:
        t_stat, p_value = stats.kruskal(*[v for v in values_distance])
        logger.info("Executive Kruskal test")
        #dict_p_value[f'{pheno_from}to{pheno_to}']=p_value
    #print(p_value)
    return p_value
    



def plot_distance_curve(path_output_result,df,pheno_from,pheno_to,p_value):
    dire =os.path.join(path_output_result,"distance_curve")
    if not os.path.exists(dire):
        os.makedirs(dire)

    #lista dei gruppi presenti
    groups=list(df["GROUP"].unique())

    #array contenente i valori delle distanze per le diverse cellulle per i diversi pazienti nei diversi gradi [[valori grado II],[valori grado III]]
    values_distance=[]

    for g in groups:
        temp=list(df[df["GROUP"]==g]["DISTANCE"].values)
        values_distance.append(temp)
    
    #group_labels = df["GROUP"].unique()
    #fig = ff.create_distplot(values_distance, group_labels,show_hist=False,show_rug=False)
    fig = ff.create_distplot(values_distance, groups,show_hist=False,show_rug=False)
    
    if p_value < 0.05 and p_value >= 0.001:
        fig.update_layout({'plot_bgcolor':'white'},title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:green;'>p value < 0.05</span>",
                        yaxis=dict(tickformat=".4f",title_text=r"$Density$"),xaxis=dict(title_text=r"$Distance_{z}$"),legend_title_text="Group")
        fig.write_image(f'{dire}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.png',scale=6)
    
    elif p_value < 0.001:
        fig.update_layout({'plot_bgcolor':'white'},title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:green;'>p value < 0.001</span>",
                        yaxis=dict(tickformat=".4f",title_text=r"$Density$",gridcolor='lightgrey'),xaxis=dict(title_text=r"$Distance_{z}$",gridcolor='lightgrey'),legend_title_text="Group")
        fig.write_image(f'{dire}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.png',scale=6)
        
    elif p_value >= 0.05 and p_value<10:
        fig.update_layout({'plot_bgcolor':'white'},title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:red;'>p value > 0.05</span>",
                        yaxis=dict(tickformat=".4f",title_text=r"$Density$",gridcolor='lightgrey'),xaxis=dict(title_text=r"$Distance_{z}$",gridcolor='lightgrey'),legend_title_text="Group")
        fig.write_image(f'{dire}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.png',scale=6)

def main():
    
    logger.info("######################## STATISTICAL ANALYSIS DISTANCE #########################")
    
    logger.info("Reading configuration file")
    with open("config.json") as f:
        data=json.load(f)
    
    #path folder with distance
    root_folder=os.path.join(data["Paths"]["output_folder"],"Distance")

    if not os.listdir(root_folder):
        logger.error(f"{root_folder} is an empty directory")
        exit()
    
    #path folder save distance statistical results
    path_output=os.path.join(data["Paths"]["output_folder"],"Distance_Statistical")

    #create Distance_Statistical folder
    create_output_dir(path_output)
   
    #list of pheno of interested
    pheno_interested=data["Phenotypes"]["pheno_list"]

    #pheno from if presents
    pheno_from=data["statistical_distance"]["pheno_from"]

    #pheno_to if presents
    pheno_to=data["statistical_distance"]["pheno_to"]

    #flag plot distance curve
    plot_distance=data["statistical_distance"]["plot_distance"]

    #flag save csv of zeta score values
    save_csv_zetascore=data["statistical_distance"]["save_csv_zetascore"]

    #groups of analysis
    groups=[]
    for group in [f for f in os.listdir(root_folder) if not f.startswith('.')]:
        groups.append(group)
    
    #path for summary statistical file
    path_stats_file=os.path.join(path_output,"summary_statistical_2.csv")

    #delete file if present
    
    if os.path.exists(path_stats_file):
        os.remove(path_stats_file)

    #open and write file in append mode
    #if len(groups)!=1:
       # f_stat=open(path_stats_file,"a")
        #f_stat.write("Phenotype"+"\t"+"Distance_to"+"\t"+"P_value"+"\t"+"\t".join(f"Meadian_{gruppo}" for gruppo in groups)+"\tGroup_major"+"\n")


    #condition if there is a pheno_from and a pheno_to of interest
    if pheno_from=="" and pheno_from=="":

        dict_statistical_result={}

        #permutation of phenotype
        for pheno in itertools.permutations(pheno_interested,2):
            pheno_from=pheno[0]
            pheno_to=pheno[1]

            logger.info(f"Analysis for distance from {pheno_from}---to---{pheno_to}")
        
            dict_distance=prepare_dataframe_distances(root_folder,pheno_from,pheno_to)

            grade_major,dict_median=calculate_median_distribution(dict_distance,groups)
           
            df_distance=create_df_distances(dict_distance,path_output,pheno_from,pheno_to,save_csv_zetascore)

            if len(df_distance)==0:
                logger.warning(f"No Distance for {pheno_from}--{pheno_to}")
                continue

            ##BUG è NECESSARIO RILEGGERLO????
            #df_distance=pd.read_csv(f"{path_output}/csv/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv",sep="\t")
            
            box_plots_distances(path_output,df_distance,pheno_from,pheno_to)
            
            if not f"{pheno_from}to{pheno_to}" in dict_statistical_result.keys():
                dict_statistical_result[f"{pheno_from}to{pheno_to}"]={}

            pvalue=statistical_test(path_output,df_distance,pheno_from,pheno_to)

            dict_statistical_result[f"{pheno_from}to{pheno_to}"]["p_value"]=pvalue
            dict_statistical_result[f"{pheno_from}to{pheno_to}"]["grade_major"]=grade_major
            dict_statistical_result[f"{pheno_from}to{pheno_to}"]["median"]=dict_median

            if plot_distance:
                plot_distance_curve(path_output,df_distance,pheno_from,pheno_to,pvalue)


        with open(path_stats_file,"w") as f:
            f.write("Phenotype"+"\t"+"Distance_to"+"\t"+"P_value"+"\t"+"\t".join(f"Meadian_{gruppo}" for gruppo in groups)+"\tGroup_major"+"\n")
            for pheno, stat in dict_statistical_result.items():
                pheno=pheno.split("to")
                pheno_from=pheno[0]
                pheno_to=pheno[1]
                p_value=stat["p_value"]
                grade_major=stat['grade_major']
                f.write(pheno_from+"\t"+pheno_to+"\t"+str(p_value)+"\t"+"\t".join(str(stat["median"][g]) for g in groups)+"\t"+grade_major+"\n")

            #statistical_curve_plot(path_output,df_distance,pheno_from,pheno_to,grade_major,dict_median,f_stat,groups,plot_distance)
    

    #if there are a pheno_from and a pheno-to of interest
    else:
        dict_statistical_result={}

        if pheno_from not in pheno_interested or pheno_to not in pheno_interested:
            logger.critical(f"Name error inserted in pheno from and/or pheno to ({pheno_from}---{pheno_to})")
            return()

        dict_distance=prepare_dataframe_distances(root_folder,pheno_from,pheno_to)
        grade_major,dict_median=calculate_median_distribution(dict_distance,groups)
        df_distance=create_df_distances(dict_distance,path_output,pheno_from,pheno_to,save_csv_zetascore)
        if len(df_distance)==0:
            logger.critical(f"No Distance for {pheno_from}--{pheno_to}")
            return()
            
        #df_distance=pd.read_csv(f"{path_output}/csv/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv",sep="\t")
        box_plots_distances(path_output,df_distance,pheno_from,pheno_to)

        if not f"{pheno_from}to{pheno_to}" in dict_statistical_result.keys():
                dict_statistical_result[f"{pheno_from}to{pheno_to}"]={}

        pvalue=statistical_test(path_output,df_distance,pheno_from,pheno_to)

        dict_statistical_result[f"{pheno_from}to{pheno_to}"]["p_value"]=pvalue
        dict_statistical_result[f"{pheno_from}to{pheno_to}"]["grade_major"]=grade_major
        dict_statistical_result[f"{pheno_from}to{pheno_to}"]["median"]=dict_median

        if plot_distance:
                plot_distance_curve(path_output,df_distance,pheno_from,pheno_to,pvalue)
        
        with open(path_stats_file,"w") as f:
            f.write("Phenotype"+"\t"+"Distance_to"+"\t"+"P_value"+"\t"+"\t".join(f"Meadian_{gruppo}" for gruppo in groups)+"\tGroup_major"+"\n")
            for pheno, stat in dict_statistical_result.items():
                pheno=pheno.split("to")
                pheno_from=pheno[0]
                pheno_to=pheno[1]
                p_value=stat["p_value"]
                grade_major=stat['grade_major']
                f.write(pheno_from+"\t"+pheno_to+"\t"+str(p_value)+"\t"+"\t".join(str(stat["median"][g]) for g in groups)+"\t"+grade_major+"\n")



    #if not os.listdir(path_output):
       # logger.info("delete empty folder Distance_Statistical")
       # os.rmdir(path_output)
    logger.info("End statistical analysis!")      
    return()
    
if __name__=="__main__":
    main()
