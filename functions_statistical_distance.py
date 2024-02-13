import logging
logging.basicConfig(format='[%(levelname)s] ALOA - %(asctime)s - %(message)s',level=logging.DEBUG)

def create_folder_dist_standard(root_folder):
    import os
    import pathlib
    
    pathlib.Path(os.path.join(root_folder,"Stats")).mkdir(parents=True,exist_ok=True)

def standardization_distance_all_image(values,paz):
    if len(values)<2:
        print(f"only one or less distance for {paz}")
        return [],0,0

    import numpy as np
    mean=np.mean(values)
    std=np.std(values)

    if std==0:
        print(f"for {paz} std = 0")
        return [],0,0
                        
    values_standard= list(map(lambda x: (x - mean)/std,values))
    return values_standard,mean,std


def standardization_distance_single_roi(values,root_folder,pheno_from,pheno_to,paz):
    import pandas as pd
    import os
    file_stats=os.path.join(root_folder,"Stats","stats.csv")
    df_stats=pd.read_csv(file_stats,sep="\t")
    #print(paz)
    #print(df_stats.loc[(df_stats["ID"]==paz) & (df_stats["Phenotype"]==pheno_from) & (df_stats["Distance_to"]==pheno_to)]["Mean"])
    if len(values)<2:
        print(f"only one distance for {paz}")
        return []
    mean=df_stats.loc[(df_stats["ID"]==paz) & (df_stats["Phenotype"]==pheno_from) & (df_stats["Distance_to"]==pheno_to)]["Mean"]
    std=df_stats.loc[(df_stats["ID"]==paz) & (df_stats["Phenotype"]==pheno_from) & (df_stats["Distance_to"]==pheno_to)]["Std"]
    print(paz,pheno_from,pheno_to,std)

    if mean.empty or float(std)==0:
        print(f"std = 0 for {paz}")
        return []
    
    values_standard= list(map(lambda x: (x - float(mean))/float(std),values))
    return values_standard

    



def prepare_dataframe_distances(root_folder,pheno_from,pheno_to,f,all_roi=True):
    import polars as pl
    import os
    import numpy as np
    import math

    DICT_DISTANCES = dict()

    list_groups = os.listdir(root_folder)
    for group in list_groups:
        if group==".DS_Store":
            continue
        if group=="Stats":
            continue
        print("Parsing Group:",group)
        DICT_DISTANCES[group] = []
        for directory in os.listdir(os.path.join(root_folder,group)):
            if directory==".DS_Store":
                continue
            if directory=="csv":
                list_files=os.listdir(os.path.join(root_folder, group,directory))
                for file in list_files:
                    paz=file.split("_")[0]
                    file_path = os.path.join(root_folder, group,directory, file)

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
                        #print("COLONNA NON TROVATA PER IL FILE", file)
                        break
                    #print(file, "NOME COLONNA:", real_pheno_to)

                    #PRENDI DATI
                    df_filtered = pl.scan_csv(
                            file_path,
                            separator="\t",
                            null_values=["NA"]
                        ).select(["Phenotype", real_pheno_to]) \
                        .drop_nulls() \
                        .filter(pl.col("Phenotype") == pheno_from) \
                        .select([real_pheno_to])
                    
                    
                    values = list(map(lambda x: float(x),df_filtered.collect()[real_pheno_to].to_list()))
                    values = [item for item in values if not(math.isnan(item)) == True]
                    if all_roi:
                        try:
                            values_standard,mean,std= standardization_distance_all_image(values,paz)
                            if values_standard !=[]:
                                f.write(f'{paz}\t{pheno_from}\t{pheno_to}\t{mean}\t{std}\n')
                            #DICT_DISTANCES[group] += values
                            DICT_DISTANCES[group] += values_standard
                        except:
                            continue
                        
                    else:
                        values_standard=standardization_distance_single_roi(values,root_folder,pheno_from,pheno_to,paz)
                        DICT_DISTANCES[group] += values_standard

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
    import os
    
    temp_folder=os.path.join(path_output_results,"Distance_Statistical")
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

        #---->logging
        logging.info(f"Created folder {temp_folder}")
        
    else:
        #--->logging
        logging.info(f"Folder {temp_folder} already exists")
        


#*****************************************************************
def create_df_distances(DICT_DISTANCES,path_output,pheno_from,pheno_to):
    import pandas as pd
    import os

    pre_df = []
    for k, v in DICT_DISTANCES.items():
        for e in v:
            pre_df.append((k, e))

    df = pd.DataFrame(pre_df, columns=["GROUP", "DISTANCE"])
    if not os.path.exists(path_output):
        os.mkdir(path_output)
        df.to_csv(path_output, sep="\t",index=False)
    if len (df)!=0:
        direc=os.path.join(path_output,"csv")
        if not os.path.exists(direc):
            os.makedirs(direc)
            df.to_csv(f'{direc}/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv', sep="\t",index=False)
        df.to_csv(f'{direc}/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv', sep="\t",index=False)
        #df.to_csv(f'{path_output}/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv', sep="\t",index=False)
    else:
        logging.warning(f"dataframe is empty for distance from {pheno_from} to {pheno_to}")

    return df


def calculate_median_distribution(dictionary_group,pheno_from,pheno_to,path_result):
    import numpy as np
    import pandas as pd
    import os
    import math
    dict_mean={}
    
    dict_mean["Grade_II"]=None
    dict_mean["Grade_III"]=None

    
    for _k,_v in dictionary_group.items():
        if len (_v)==0:
            continue
        mean=np.median(_v)
        #if math.isnan(mean):
            #continue
            #if _k not in dict_mean.keys():
                #dict_mean[_k]=0
        dict_mean[_k]=mean
        
    grade_major=""
    if dict_mean["Grade_II"] is None or dict_mean["Grade_III"] is None:
        print("for one o both grade value is none")
        return 
    elif dict_mean["Grade_II"]>dict_mean["Grade_III"]:
        grade_major="Grade_II"
    else:
        grade_major="Grade_III"

    dire=os.path.join(path_result,"Comparison_distance_mean")
    if not os.path.exists(dire):
        os.makedirs(dire)
            
    with open(os.path.join(dire,f"distance_mean_{pheno_from}_to_{pheno_to}.csv"),"w") as f:
        f.write("Mean_Grade_II\tMean_Grade_III\tResults\n")
        f.write(str(dict_mean["Grade_II"])+"\t"+str(dict_mean["Grade_III"])+"\t"+grade_major+"\n")

    



    




#*****************************************************************

def box_plots_distances(path_ouput_results,df,pheno_from,pheno_to):
    import tap
    import os
    dire=os.path.join(path_ouput_results,"box_plot")
    if not os.path.exists(dire):
        os.makedirs(dire)
    if len(df["GROUP"].unique())!=1: 
        tap.plot_stats(df,x="GROUP",y="DISTANCE",filename=f'{dire}/distances_box_plot_{pheno_from}_to_{pheno_to}.png',kwargs={"title":f'Distance from {pheno_from} to {pheno_to}',"labels":{"GROUP": "Group",
                     "DISTANCE": r'$Distance_{z}$',
                     "GROUP": 'Group'
                 }})

#*****************************************************************

def statistical_curve_plot(path_output_result,df,pheno_from,pheno_to):
    from scipy import stats
    import seaborn as sns
    import matplotlib.pyplot as plt
    import plotly.figure_factory as ff
    import os

    dire =os.path.join(path_output_result,"distance_curve")
    if not os.path.exists(dire):
        os.makedirs(dire)


    groups=list(df["GROUP"].unique())
    values_distance=[]

    for g in groups:
        temp=list(df[df["GROUP"]==g]["DISTANCE"].values)
        values_distance.append(temp)
    
    

    p_value=10

    if len(df["GROUP"].unique())==1:
        print("only one group")
   
        hist_data=[list(df["DISTANCE"].values)]
        
        group_labels = df["GROUP"].unique()

        fig = ff.create_distplot(hist_data, group_labels,show_hist=False,show_rug=False)
        # Add title
        fig.update_layout(title_text=f"Distance-Z score from {pheno_from} to {pheno_to}",
                          yaxis=dict(tickformat=".4f",title_text=r"$Density$"),xaxis=dict(title_text=r"$Distance_{z}$"),legend_title_text="Group")
        fig.write_image(f'{dire}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.jpeg',scale=6)
        #fig.show()

        #ax=sns.kdeplot(data=df,x="DISTANCE")
        #plt.title(f'Distance {pheno_from} to\n{pheno_to}')
        #plt.legend(df["GROUP"].values)
        #plt.savefig(f'{path_output_result}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.png')
        #plt.close()
    
    elif len(df["GROUP"].unique())==2:
        t_stat, p_value = stats.mannwhitneyu(values_distance[0], values_distance[1])

    elif len(df["GROUP"].unique())>2:#se i gruppi sono più di due, bisognerà utilizzare il test di Kruskal per più sample
        t_stat, p_value = stats.kruskal([v for v in values_distance])
    #p_value= "{:.3e}".format(p_value)
    

    #hist_data=[list(df["DISTANCE"].values)]
        
    group_labels = df["GROUP"].unique()

    fig = ff.create_distplot(values_distance, group_labels,show_hist=False,show_rug=False)
    
    
    #ax=sns.kdeplot(data=df,x="DISTANCE",hue="GROUP")
    #x=(max(df["DISTANCE"].values)/2)
    #y=(max(ax.lines[0].get_ydata())/2)
    if p_value < 0.05 and p_value >= 0.001:
        fig.update_layout(title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:green;'>p value < 0.05</span>",
                          yaxis=dict(tickformat=".4f",title_text=r"$Density$"),xaxis=dict(title_text=r"$Distance_{z}$"),legend_title_text="Group")
        fig.write_image(f'{dire}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.jpeg',scale=6)
    
        #ax.text(x, y,f"p-value < 0.05")
    elif p_value < 0.001:
        fig.update_layout(title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:green;'>p value < 0.001</span>",
                          yaxis=dict(tickformat=".4f",title_text=r"$Density$"),xaxis=dict(title_text=r"$Distance_{z}$"),legend_title_text="Group")
        fig.write_image(f'{dire}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.jpeg',scale=6)
        
    elif p_value >= 0.05 and p_value<10:
        fig.update_layout(title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:red;'>p value > 0.05</span>",
                          yaxis=dict(tickformat=".4f",title_text=r"$Density$"),xaxis=dict(title_text=r"$Distance_{z}$"),legend_title_text="Group")
        fig.write_image(f'{dire}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.jpeg',scale=6)
    
        
        
    #plt.title(f'Distance {pheno_from} to\n{pheno_to}')
    #plt.savefig(f'{path_output_result}/plot_statistical_distance_{pheno_from}_to_{pheno_to}.png')
    #plt.close()
    #fig.show()
    #plt.show()




def main():
    import json
    import itertools
    import os
    import pandas as pd

    with open("config.json") as f:
        data=json.load(f)
    

    root_folder=data["statistical_distance"]["root_folder"]
    #dir_output="./output/output_PanelA/"
    dir_output=data["Paths"]["output_folder"]
    dit_distance_stats=create_folder_dist_standard(root_folder)
    #path_output="./output/output_PanelA/Distance_Statistical/"
    path_output=data["statistical_distance"]["save_folder"]
    #pheno_interested=list(data["Descriptive"]["pheno_interested"].keys())
    pheno_interested=data["Phenotypes"]["pheno_list"]
    all_roi=data["statistical_distance"]["all_roi"]
    #print(pheno_interested)
    pheno_from=data["statistical_distance"]["pheno_from"]
    pheno_to=data["statistical_distance"]["pheno_to"]

    with open(os.path.join(root_folder,"Stats","stats.csv"),"a") as f_stat:
        f_stat.write("ID\tPhenotype\tDistance_to\tMean\tStd\n")
        if pheno_from=="" and pheno_from=="":
            for pheno in itertools.permutations(pheno_interested,2):
                pheno_from=pheno[0]
                pheno_to=pheno[1]
                logging.info(f"{pheno_from}---to---{pheno_to}")
                #print(pheno_from,"---to---",pheno_to)
            
                dict_distance=prepare_dataframe_distances(root_folder,pheno_from,pheno_to,f_stat,all_roi)
                calculate_median_distribution(dict_distance,pheno_from,pheno_to,path_output)

                create_output_dir(dir_output)
                df_distance=create_df_distances(dict_distance,path_output,pheno_from,pheno_to)
                if len(df_distance)==0:
                    logging.warning(f"No Distance for {pheno_from}--{pheno_to}")
                    #print(f"----WARNING----\nNo Distance for {pheno_from}--{pheno_to}")
                    dire=os.path.join(path_output,"no_distance")
                    if not os.path.exists(dire):
                        os.makedirs(dire)
                    with open (f"{dire}/no_distance_matching.txt","a") as f:
                        f.write(f"from {pheno_from} --- to {pheno_to}\n")
                    continue
                df_distance=pd.read_csv(f"{path_output}/csv/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv",sep="\t")
                box_plots_distances(path_output,df_distance,pheno_from,pheno_to)
                #print("si è rotto", df_distance.describe())
                statistical_curve_plot(path_output,df_distance,pheno_from,pheno_to)
        
        
        else:
            dict_distance=prepare_dataframe_distances(root_folder,pheno_from,pheno_to,f_stat,all_roi)
            calculate_median_distribution(dict_distance,pheno_from,pheno_to,path_output)
            create_output_dir(dir_output)
            df_distance=create_df_distances(dict_distance,path_output,pheno_from,pheno_to)
            if len(df_distance)==0:
                logging.warning(f"No Distance for {pheno_from}--{pheno_to}")
                #print(f"----WARNING----\nNo Distance for {pheno_from}--{pheno_to}")
                dire=os.path.join(path_output,"no_distance")
                if not os.path.exists(dire):
                    os.makedirs(dire)
                with open (f"{dire}/no_distance_matching.txt","a") as f:
                    f.write(f"from {pheno_from} --- to {pheno_to}\n")
                #continue
                
            df_distance=pd.read_csv(f"{path_output}/csv/df_statistical_distance_{pheno_from}_to_{pheno_to}.csv",sep="\t")
            box_plots_distances(path_output,df_distance,pheno_from,pheno_to)
            #print("si è rotto", df_distance.describe())
            statistical_curve_plot(path_output,df_distance,pheno_from,pheno_to)

            

if __name__=="__main__":
    main()
