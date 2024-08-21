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

def raw_count_cells(PATH_MERGE_FOLDER,list_pheno):
    '''
    function to calculate raw counts from merged file
    Parameters
    ----
    PATH_MERGE_FOLDER : str
    list_pheno : array

    Return
    ----
    dict_global_info : dict
    '''

    #creating a dictionary with raw count information for each patient
    
    logger.info("Start raw count calculation:")

    dict_global_info = {}

    #iterating on folders presented in path merge folder
    for _d in [f for f in os.listdir(PATH_MERGE_FOLDER) if not f.startswith('.')]:
        sub_directory=os.path.join(PATH_MERGE_FOLDER,_d)
        if not os.listdir(sub_directory):
            logger.error(f"Directory {_d} is empty")
            exit()

        #iterating on each files presented into group folder
        for f in [f for f in os.listdir(sub_directory) if not f.startswith('.')]:
            
            #patient name
            _id_paz=f.replace("Merge_cell_seg_data_clean_","").replace(".txt","")
            logger.info(f"Subject: {_id_paz}")
            #adding group name as key of dictionary
            if _d not in dict_global_info.keys():
                dict_global_info[_d]={}

            #adding patient name as second key of dictionary
            if _id_paz not in dict_global_info[_d].keys():
                dict_global_info[_d][_id_paz]={}

            #path of patient merged file
            _file=os.path.join(sub_directory,f)
    
            #reading the file
            _data=pd.read_csv(_file, sep="\t")
            _data["Pheno"]=list(map(lambda x: x.replace("+","+,"),_data["Pheno"]))
            _data["Pheno"]=list(map(lambda x: x[:-1],_data["Pheno"]))

            #---->logging
            logger.info(f"Reading File {_file}")

            #calculating total cells
            _total_cells=len(_data)

            #adding Total_Cells as third key of dictionary
            dict_global_info[_d][_id_paz]["Total_Cells"]=_total_cells

            #count of pheno column
            _data_groupped = _data.groupby(["Pheno"])["Pheno"].count()

            #eliminating other from pheno column
            #_data_dict = {frozenset(k.replace("other,", "").replace(",other", "").split(",")): v for k, v in _data_groupped.to_dict().items()}
            _data_dict = {frozenset(k.replace("OTHER,", "").replace(",OTHER", "").split(",")): v for k, v in _data_groupped.to_dict().items()}
            
            #adding phenotype as fourth key of dictionary
            for main_pheno in list_pheno:
                dict_global_info[_d][_id_paz][f"{main_pheno}"] = _data_dict.get(frozenset(main_pheno.split(",")), 0)

    logger.info("End raw count calculation")

    return dict_global_info 

#******************************************************

def calculate_mean_group_cells(dictionary_raw_count):
    '''
    function to calculate the mean of the total cells for each groups

    Parameters
    ----
    dictionary_raw_count: dict

    Return
    ----
    count_mean_group_cells: dict
    '''
    
    #creating dictionary with the mean of total cells for each group
    count_mean_group_cells={}

    for _k,_v in dictionary_raw_count.items():
        _temp=[]
        #adding group name as key of the dictionary
        if _k not in count_mean_group_cells.keys():
            count_mean_group_cells[_k]=0
        for _,_v2 in _v.items():
            #adding on _temp array the total cells of each patient for each group
            _temp.append(_v2["Total_Cells"])
        
        #calculating the mean of each group
        mean_group=round(np.mean(_temp),2)

        #adding the mean as value of the dictionary
        count_mean_group_cells[_k]=mean_group

    return count_mean_group_cells

#******************************************************

def normalized_count_cells(dictionary_raw_count,dictionary_mean_count):
    '''
    - formula normalizzation single group = (number of cells positive for target (or phenotipes) / total cells patient ) * mean of total cells of group

    function to calculate normal count for each patient

    Parameters
    ----
    dictionary_raw_count : dict
    dictionary_mean_count : dict

    Return
    ----
    dictionary_normalized_count_cells : dict

    '''

    logger.info("Start normalized count calculation")

    #creating dictionary for normalized count
    dict_normalized_count_cells={}

    #iterating on raw count dictionary
    for group, patients in dictionary_raw_count.items():
        # adding group as first key of the dictionary
        dict_normalized_count_cells[group]={}

        # taking the mean value of the group
        mean_group=dictionary_mean_count[group]

        #iterationg on patient- pheno raw count dictionary
        for patient, counters in patients.items():
            logger.info(f"Subject: {patient}")

            #adding patient as second key of the dictionary
            dict_normalized_count_cells[group][patient]={}

            #taking the value of total cells for each patient
            total_cells=counters["Total_Cells"]
          
            #adding Total_Cells as third key of the dictionary
            dict_normalized_count_cells[group][patient]["Total_Cells"]=total_cells

            #iterating on pheno-raw count
            for pheno, value in counters.items():
                if pheno == "Total_Cells":
                    continue

                # calculating the normalized count 
                _new_value=(value/total_cells)*mean_group

                #adding pheno as fourth key of the dictionary
                dict_normalized_count_cells[group][patient][pheno]=_new_value


    logger.info("End normalized count calculation")

    return dict_normalized_count_cells

#******************************************************

def create_output_dir(path_output,groups_name):
    '''
    function to create the folder where descriptive's output are saved

    Parameters
    ----
    path_output : str
    groups_name : list

    Return
    ----
    None

    '''
    
    for _g in groups_name:
        temp_folder=os.path.join(path_output,_g)
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
            logger.info(f" Created folder {temp_folder}")
        
#******************************************************

def create_summury_file(path_output_results,dictionary_count,type_data):
    '''
    Function to create a csv file for each patient, for each group, where are saved the raw count and/or normalized count for the phenotypes of interest
    
    Parameters
    ----
    path_output_results : str
    dictionary_count : dict
    type_data : str

    Return
    ----
    None
    '''
    logger.info("writing csv file with count info for each patient")
    #iterationg on dictionary count
    for group,patients in dictionary_count.items():
        logger.info(f"Group: {group}")
        #path of csv directory into group directory
        dire = os.path.join(path_output_results,group,"csv")
        try:
            pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        except Exception:
            logger.critical(f"Something went wrong during the creation of {dire} folder!")
            return()
        #iterating on patient and pheno count
        for patient,phenos in patients.items():
            logger.info(f"Subject: {patient}")

            with open(f'{dire}/{type_data}_count_{patient}.csv',"w") as f:
                f.write(f"Patient\tPheno\tCount_{type_data}\n")

                #----> aggiunta del log per l'inizio scrittura del file
                logger.info(f"Start Writing csv file...")

                #list of phenotypes
                pheno=list(phenos.keys())
                #iteration on pheno list
                for p in pheno:
                    if p=="Total_Cells":
                        f.write(f'{patient}\t{p}\t{dictionary_count[group][patient]["Total_Cells"]}\n')
                    else:
                        f.write(f'{patient}\t{p}\t{dictionary_count[group][patient][p]}\n')
                
                #----> aggiunta del log per la fine della scrittura del file
                logger.info("End Writing " + os.path.join(dire,type_data + "_count_" + patient + ".csv") + " file for " + patient + " patient")

#***********************************+

def create_norm_all_file(path_output_results,dictionary_norm_count):

    '''
    
    Function to save the normalized count, for all groups that are present, into a csv file

    Parameters
    ----
    path_output_results : str
    dictionary_norm_count : dict

    Return
    ----
    None

    '''

    logger.info("Start normalized count calculatin on all groups")

    for group,patients in dictionary_norm_count.items():
        dire = os.path.join(path_output_results,group,"csv")
        try:
            pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        except Exception:
            logger.critical(f"Something went wrong during the creation of {dire} folder!")
            return()

        for patient,phenos in patients.items():
            with open(os.path.join(dire, "all_norm_count_"+patient+".csv"),'w') as f:
                f.write(f"Patient\tPheno\tCount_all_Norm\n")

                logger.info(f"Writing normalized count file on all groups for patient {patient}")
                pheno=list(phenos.keys())
                for p in pheno:
                    if p=="Total_Cells":
                        f.write(f'{patient}\t{p}\t{dictionary_norm_count[group][patient]["Total_Cells"]}\n')
                    else:
                        f.write(f'{patient}\t{p}\t{dictionary_norm_count[group][patient][p]}\n')

    logger.info("End normalized count calculation on all groups")

#******************************************************

def bar_plot(path_output_result,dict_data,type_data):
    '''
    Function to create bar plot, for each group, with raw and/or normalized count of single target or phenotype of interest

    Parameters
    ----
    PATH_OUTPUT_RESULTS : str
    dictionary_count : dict
    type_data : str

    Return
    ----
    None

    '''
    
    logger.info(f"Barplot creation for {type_data} count")
    
    for group,data in dict_data.items():
        
        dire = os.path.join(path_output_result,group,"Bar_plot")

        #check if the directory exists
        try:
            pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        except Exception:
            logger.critical(f"Something went wrong during the creation of {dire} folder!")
            return()
        
        logger.info(f"Group: {group}")

        max_val_y=0
        fig=go.Figure()
        for patient,targets in data.items():

            logger.info(f"Subject: {patient}")

            #sorted phenotypes in alphabetic order
            sorted_target=dict(sorted(targets.items(),key=lambda x:x[0]))
            
            #transforming pheno into a list
            target=list(sorted_target.keys())

            #removing Total_Cells as pheno and its value
            target.remove("Total_Cells")
            values=list(sorted_target.values())[:-1]

            #define max value for y axes limit
            max_val=max(values)
            if max_val>max_val_y:
                max_val_y=max_val
            # create figure
            fig.add_trace(go.Bar(x=target,y=values, name=patient))

        #adding layout on figure
        fig.update_layout(barmode='group',title_text=f"Phenotypes {type_data}-{group}",xaxis=dict(
        title="Phenotypes",showticklabels=True),yaxis=dict(title=f'log({type_data} Count)',showticklabels=True,type="log"),legend_title_text="Patients")

        #adding y axis range in log scale
        fig.update_yaxes(range=[0,(math.log(max_val_y,10)+1.0)])

        # saved figure in jpeg format
        fig.write_image(f'{dire}/Bar_Plot_{type_data}.jpeg',scale=6)
        
#******************************************************

def calculate_mean_total_groups(dictionary_raw_count):
    '''
    function to calculate the mean of the total cells, considering all groups

    Parameters
    ----
    dictionary_raw_count: dict

    Return
    ----
    mean_cells : float
    '''
  
    _count_total=[]

    #iterating on raw count dictionary
    for _k,_v in dictionary_raw_count.items():
        for _k2,_v2 in _v.items():
            # adding to _count_total list the values of total cells of each patient of all groups
            _count_total.append(_v2["Total_Cells"])

    #calculatin the global mean
    mean_cells=np.mean(_count_total)
        
    return mean_cells

#******************************************************

def normalized_count_on_all_groups(dictionary_raw_count,mean_all_groups):
    '''
    formula normalizzation all groups = number of cells positive for target (or phenotipes) / total cells patient ) * mean of total cells of all groups

    Functions to calculate the normalized count on all groups present

    Parameters:
    ----
    dictionary_raw_count : dict
    mean_cells : float

    Return:
    ----
    dict_normalized_count_cells : dict

    '''

    logger.info("Starting normalized count on all groups")

    # define dictionary of count normalized on all group present
    dict_normalized_count_cells={}

    #iterating on raw count dictionary
    for group, patients in dictionary_raw_count.items():
        logger.info(f"Group: {group}")
        # adding group as first key of dictionary
        dict_normalized_count_cells[group]={}
        # iterating on patient-raw count
        for patient, counters in patients.items():
            logger.info(f"Subject: {patient}")
            # adding patient as second key of dictionary
            dict_normalized_count_cells[group][patient]={}
            # save total cells of each patient
            total_cells=counters["Total_Cells"]
            # adding Total_Cells as third key of dictionary
            dict_normalized_count_cells[group][patient]["Total_Cells"]=total_cells
            # iterating on phenotypes and relative count
            for pheno, value in counters.items():
                if pheno == "Total_Cells":
                    continue
                #normalized count base on mean of all groups
                _new_value=(value/total_cells)*mean_all_groups
                # adding phenotypes as fourth key of dictionary
                dict_normalized_count_cells[group][patient][pheno]=_new_value

    logger.info("End normalized count on all groups")

    return dict_normalized_count_cells

#******************************************************

def prepare_data_box_plot(dictionary_count):
    '''
    Function to prepare data for next analysis, when there are of 2 or more groups. Creates a Dataframe with Group, Target, Value columns.

    Parameters:
    ----
    PATH_MERGE_FOLDER : str
    dictionary_count : dict

    Return
    ----
    _data_all: DataFrame

    '''

    #creating an empty dataframe
    _data_all=pd.DataFrame()
    _data_norm=[]

    # iterating on dictionary count
    for _k,_v in dictionary_count.items():
        for _,_v2 in _v.items():
            for _t in _v2.keys():
                if _t=="Total_Cells":
                    continue

                #array with information about group, phenotype and count
                data_temp=[_k,_t,_v2[_t]]
                _data_norm.append(data_temp)

    #append list on dataframe
    _data_all=pd.DataFrame(_data_norm,columns=["group","pheno","Count"])

    return _data_all
  
#*****************************************************************

def create_output_box_plot_dir(path_output):
    '''
    Function to create the Box Plots folder, where are saved box plot image
    Funzione per la creazione cartella dove veranno inseriti gli output, con creazione della cartella box plot dove saranno salvate le immagini

    Parameters
    ----
    path_outpu : str

    Return
    ----
    None

    '''
    
    temp_folder=os.path.join(path_output,"Box Plots")
    
    try:
        pathlib.Path(temp_folder).mkdir(parents=True, exist_ok=True)
    except Exception:
        logger.critical(f"Something went wrong during the creation of {temp_folder} folder!")
        return()
    
#*****************************************************************

def create_comparison_box_plot(path_output_result,data_all,p_adjust,test,type_data):
    #print(data_all)
    '''
    Function to make statistical test and annoteted the relative boxplots, generated with ploylt, to compare  phenotypes betweenness different group

    Parameters
    ----
    path_output : str
    data_all : DataFrame
    type_data : str
    

    Return
    ---
    None
   
    '''

    logger.info(f"Box Plot creation on {type_data} count")

    #define paramerts as x, y and separated
    x="pheno"
    y="Count"
    hue="group"
    hue_order=list(data_all["group"].unique())
    lables=sorted(list(data_all["pheno"].unique()))

    filename=os.path.join(path_output_result,"Box Plots","box_plot_comparison_"+type_data+".jpeg")

    if len(hue_order)==2 and test!="paired":
        logger.info(f"Found {len(hue_order)} groups: Start Mann-Whitney Test")
        try:
            tap.plot_stats(data_all,x,y,order=lables,filename=filename,type_correction=p_adjust,export_size=(1400, 950, 3),subcategory=hue,kwargs={"width":4000,"height":1000,"title":f"Phenotypes Comparison between {hue_order[0]} and {hue_order[1]}-{type_data} Count","log_y":True,"labels":{"pheno":"Phenotypes","value":f"log({type_data} Counts)","group":"Group"}})
            logger.info("Mann-Whitney Test done successfully!")        
        except ValueError:
            logger.error("Mann-Whitney Test Error-All numbers are identical")
            return()
        except Exception as e:
            logger.error("Something went wrong during Mann-Whitney test")
            return()
    if len(hue_order)==2 and test=="paired":
        logger.info(f"Found {len(hue_order)} groups: Start Wilcoxon Test")
        try:
            tap.plot_stats(data_all,x,y,order=lables,filename=filename,type_test="wilcoxon",type_correction=p_adjust,export_size=(1400, 950, 3),subcategory=hue,kwargs={"width":4000,"height":1000,"title":f"Phenotypes Comparison between {hue_order[0]} and {hue_order[1]}-{type_data} Count","log_y":True,"labels":{"pheno":"Phenotypes","value":f"log({type_data} Counts)","group":"Group"}})
            logger.info("Wilcoxon Test done successfully!")        
        except ValueError:
            logger.error("Wilcoxon Test Error-All numbers are identical")
            return()
        except Exception as e:
            logger.error("Something went wrong during Wilcoxon test")
            return()
        
    elif len(hue_order)>2:
        logger.info(f"Found {len(hue_order)} groups: Start Kruskal-Wallis Test")
        try:
            if p_adjust is None:
                p_adjust="bonferroni"
            tap.plot_stats(data_all,x,y,type_test="dunn",order=lables,subcategory=hue,filename=filename,type_correction=p_adjust,export_size=(1400, 950, 3),kwargs={"width":4000,"height":1000,"title":f"Phenotypes Comparison between {','.join(hue_order)}-{type_data} Count","log_y":True,"labels":{"pheno":"Phenotypes","value":f"log({type_data} Counts)","group":"Group"}})
            logger.info("Kruskal-Wallis Test done successfully!")
        except ValueError:
            logger.error("Kruskal-Wallis Test Error-All numbers are identical")
            return()
        except Exception:
            logger.error("Something went wrong during Kruskal-Wallis test")
            return()

#*****************************************************************
#*****************************************************************
def main(data):
    
    print("\n############################# DESCRIPTIVE ANALYSIS #############################\n")
    
    logger.info("Start overview process: This step will give a description of input data and, if more than one group is involved, will also provide a statistical comparison for all the listed phenotypes\n")
    
    logger.info("Reading configuration file")
    
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
    
    p_adjust=data["Stats"]["p_adj"]
    p_adjust=p_adjust.lower()
    if p_adjust=="":
        p_adjust=None
    
    test=data["Stats"]["sample_type"]

    if data["Descriptive"]["raw"]:
        create_summury_file(path_output_results,dict_raw_count,"Raw")
        bar_plot(path_output_results,dict_raw_count,"Raw")
        if len(groups)>=2:
            create_output_box_plot_dir(path_output_results)
            df_raw=prepare_data_box_plot(dict_raw_count)
            create_comparison_box_plot(path_output_results,df_raw,p_adjust,test,"Raw")
        else:
            logger.critical("Found only one group - impossible to continue with groups comparison and box plot figure")
            #return()
            
    if data["Descriptive"]["normalized"]:
        mean_group=calculate_mean_group_cells(dict_raw_count)
        print(mean_group)
        norm_count=normalized_count_cells(dict_raw_count,mean_group)
        create_summury_file(path_output_results,norm_count,"Norm")
        bar_plot(path_output_results,norm_count,"Normalized")
        if len(groups)>=2:
            mean_group_all = calculate_mean_total_groups(dict_raw_count)
            norm_count_all_groups = normalized_count_on_all_groups(dict_raw_count,mean_group_all)
            create_norm_all_file(path_output_results,norm_count_all_groups)
            create_output_box_plot_dir(path_output_results)
            df_norm=prepare_data_box_plot(norm_count_all_groups)
            create_comparison_box_plot(path_output_results,df_norm,p_adjust,test,"Normalized")
        else:
            logger.critical("Found only one group-impossible to continue with group comparison")
            #return()

    logger.info("End descriptive analysis!\n")
    return()

if __name__=="__main__":
    main()