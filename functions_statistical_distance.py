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
import os
import numpy as np
import polars as pl
import math
import pandas as pd
from scipy import stats
import itertools
import plotly.figure_factory as ff
import plotly.express as px
import tap
import scikit_posthocs as sp
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def process_pheno_pair(pheno, data, root_folder, path_output, groups, plot_distance, st_test, p_adjust, save_csv_zetascore):
    pheno_from, pheno_to = pheno
    logger.info(f"Analysis for distance from {pheno_from} to {pheno_to}")
    pvalues_dict = {}  # dizionario per annotazioni barplot
    try:
        dict_distance = prepare_dataframe_distances(root_folder, pheno_from, pheno_to)
        grade_major, dict_median = calculate_median_distribution(dict_distance, groups)
        df_distance = create_df_distances(dict_distance, path_output, pheno_from, pheno_to, save_csv_zetascore)
        if len(df_distance) == 0:
            return (pheno_from, pheno_to, None, grade_major, dict_median, pvalues_dict)

        # Statistica
        pvalue, kruskal, test = statistical_test(df_distance, path_output, st_test, p_adjust)

        # Popola dizionario pvalues per i singoli fenotipi
        pvalues_dict[pheno_from] = pvalue

        


        # Barplot per confronto specifico
        df_median_group = df_distance.groupby("GROUP")["DISTANCE"].median().reset_index()
        fig_bar = px.bar(
            df_median_group,
            x="GROUP", y="DISTANCE", color="GROUP",
            title=f"Median distance barplot: {pheno_from} → {pheno_to}"
        )
        bar_dir = os.path.join(path_output, "BarPlots_all_comparisons")
        os.makedirs(bar_dir, exist_ok=True)
        fig_bar.write_image(os.path.join(bar_dir, f"barplot_{pheno_from}_to_{pheno_to}.png"))

        pairwise_pvals = box_plots_distances(path_output, df_distance, pheno_from, pheno_to, kruskal, p_adjust, test)
        save_pairwise_csv(df_distance, pheno_from, pheno_to, path_output, pairwise_pvals)


        # Curva distanza
        if plot_distance:
            plot_distance_curve(path_output, df_distance, pheno_from, pheno_to, pvalue)

        return (pheno_from, pheno_to, pvalue, grade_major, dict_median, pvalues_dict)

    except Exception as e:
        logger.error(f"Error in phenotype pair {pheno_from}-{pheno_to}: {e}")
        return (pheno_from, pheno_to, None, None, None, pvalues_dict)

def save_pairwise_csv(df_distance, pheno_from, pheno_to, path_output, pairwise_pvals):
    out_dir = os.path.join(path_output, "pairwise_pvalues")
    os.makedirs(out_dir, exist_ok=True)
    outfile = os.path.join(out_dir, f"pairwise_pvalues_{pheno_from}_to_{pheno_to}.csv")

    rows = []
    groups = sorted(df_distance["GROUP"].unique()) if not df_distance.empty else []

    # Salva solo coppie uniche (ordina le coppie)
    for i, g1 in enumerate(groups):
        for j, g2 in enumerate(groups):
            if i >= j:
                continue  # evita duplicati e diagonale
            median1 = df_distance[df_distance["GROUP"] == g1]["DISTANCE"].median() if g1 in df_distance["GROUP"].values else np.nan
            median2 = df_distance[df_distance["GROUP"] == g2]["DISTANCE"].median() if g2 in df_distance["GROUP"].values else np.nan
            pval = pairwise_pvals.get((g1, g2), np.nan)
            rows.append({
                "PHENO_FROM": pheno_from,
                "PHENO_TO": pheno_to,
                "GROUP_1": g1,
                "GROUP_2": g2,
                "MEDIAN_1": median1,
                "MEDIAN_2": median2,
                "P_VALUE": pval
            })

    df_out = pd.DataFrame(rows, columns=["PHENO_FROM", "PHENO_TO", "GROUP_1", "GROUP_2", "MEDIAN_1", "MEDIAN_2", "P_VALUE"])
    df_out.to_csv(outfile, index=False)



def plot_heatmap_and_bar(df_all, barplot_dir, phenotypes, pvalues_dict=None):
    """
    Crea heatmap NxN delle mediane delle distanze tra fenotipi,
    separate per gruppo, così da confrontare le stesse distanze tra gruppi.
    """
    # Assicurati che la cartella di output esista
    os.makedirs(barplot_dir, exist_ok=True)

    groups = df_all["GROUP"].unique()
    n_groups = len(groups)

    # Determina numero di righe/colonne per i subplot (approssimativamente quadrato)
    n_cols = min(3, n_groups)  # massimo 3 colonne per figura
    n_rows = (n_groups + n_cols - 1) // n_cols

    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=[f"Group: {g}" for g in groups]
    )

    row, col = 1, 1
    for i, group in enumerate(groups):
        df_group = df_all[df_all["GROUP"] == group]
        df_pivot = df_group.pivot_table(
            index="PHENO", columns="PHENO_TO", values="DISTANCE", aggfunc="median"
        )

        heatmap = go.Heatmap(
            z=df_pivot.values,
            x=df_pivot.columns,
            y=df_pivot.index,
            colorscale="RdBu",
            colorbar=dict(title="Median Distance"),
        )

        fig.add_trace(heatmap, row=row, col=col)

        col += 1
        if col > n_cols:
            col = 1
            row += 1

    fig.update_layout(
        height=300 * n_rows, width=400 * n_cols,
        title_text="Median distance heatmap (NxN phenotypes) per group",
    )

    fig.write_image(os.path.join(barplot_dir, "heatmap_NxN_by_group.png"))


def merge_all_pairwise_csv(path_output):
    """
    Unisce tutti i file pairwise generati in un unico CSV globale.
    """
    import glob

    folder = os.path.join(path_output, "pairwise_pvalues")
    os.makedirs(folder, exist_ok=True)

    all_files = glob.glob(os.path.join(folder, "pairwise_pvalues_*.csv"))

    if len(all_files) == 0:
        print("Nessun file pairwise trovato.")
        return

    dfs = []
    for f in all_files:
        df = pd.read_csv(f)
        dfs.append(df)

    df_final = pd.concat(dfs, ignore_index=True)

    # ordina per PHENO_FROM, PHENO_TO, GROUP_1, GROUP_2
    df_final = df_final.sort_values(
        by=["PHENO_FROM", "PHENO_TO", "GROUP_1", "GROUP_2"]
    )

    outfile = os.path.join(folder, "ALL_PAIRWISE_RESULTS.csv")
    df_final.to_csv(outfile, index=False)

    print(f"File combinato creato: {outfile}")


def standardization_distance_all_image(values,paz):
    if len(values)<2:
        logger.info(f"only one or less distance for {paz}")
        return [],0,0

    mean=np.mean(values)
    std=np.std(values)

    if std==0:
        logger.info(f"for {paz} std = 0")
        return [],0,0
                        
    values_standard= list(map(lambda x: (x - mean)/std,values))
    return values_standard,mean,std

#*****************************************************************

def prepare_dataframe_distances(root_folder,pheno_from,pheno_to):
    DICT_DISTANCES = dict()

    for group in [f for f in os.listdir(root_folder) if not f.startswith('.')]:
        if group=="Stats":
            continue
        
        DICT_DISTANCES[group] = []
        if not os.listdir(os.path.join(root_folder,group)):
            logger.error(f"{os.path.join(root_folder,group)} is an empty folder")
            exit()
        for file in [f for f in os.listdir(os.path.join(root_folder,group)) if not f.startswith('.')]:
            paz=file.replace("_Distance.txt","")
            file_path = os.path.join(root_folder, group, file)

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

            df_filtered = pl.scan_csv(
                    file_path,
                    separator="\t",
                    null_values=["NA"]
                ).select(["Phenotype", real_pheno_to]) \
                .drop_nulls() \
                .filter(pl.col("Phenotype") == real_pheno_from) \
                .select([real_pheno_to])
            
            raw_values = df_filtered.collect()[real_pheno_to].to_list()
            values = []
            for v in raw_values:
                try:
                    val = float(v)
                    if not math.isnan(val):
                        values.append(val)
                except (ValueError, TypeError):
                    logger.warning(f"Valore non numerico '{v}' trovato in {paz}, colonna {real_pheno_to} → scartato")
                    continue

            try:
                values_standard,_,_= standardization_distance_all_image(values,paz)
                DICT_DISTANCES[group] += [(paz, val) for val in values_standard]
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

        logger.info(f"Created folder {temp_folder}")
        
    else:
       logger.info(f"Folder {temp_folder} already exists-deleting its contents for the new analysis")
       for direct in os.listdir(temp_folder):
            direct_path = os.path.join(temp_folder, direct)
            if not os.path.isfile(direct_path):
                for file in os.listdir(direct_path):
                    file_path=os.path.join(direct_path,file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                os.rmdir(direct_path)

#*****************************************************************
def create_df_distances(DICT_DISTANCES, path_output, pheno_from, pheno_to, save_zetascore=False):
    pre_df = []
    for group, values in DICT_DISTANCES.items():
        for paz, val in values:
            pre_df.append((group, paz, val))

    df = pd.DataFrame(pre_df, columns=["GROUP", "PAZ", "DISTANCE"])

    # 📂 nuova sotto-cartella "Dataframes"
    df_dir = os.path.join(path_output, "Dataframes")
    os.makedirs(df_dir, exist_ok=True)

    # 1. salva tutte le cellule
    full_csv = os.path.join(df_dir, f"df_all_samples_{pheno_from}_to_{pheno_to}.csv")
    df.to_csv(full_csv, sep="\t", index=False)

    # 2. salva mediana per paziente
    if len(df) > 0:
        df_patient = df.groupby(["GROUP", "PAZ"], as_index=False)["DISTANCE"].median()
        patient_csv = os.path.join(df_dir, f"df_per_patient_{pheno_from}_to_{pheno_to}.csv")
        df_patient.to_csv(patient_csv, sep="\t", index=False)

    # 3. opzionale: salvataggio storico
    if len(df) != 0 and save_zetascore:
        dire = os.path.join(df_dir, "csv")
        os.makedirs(dire, exist_ok=True)
        df.to_csv(os.path.join(dire, f"df_statistical_distance_{pheno_from}_to_{pheno_to}.csv"), sep="\t", index=False)
    else:
        logger.warning(f"Dataframe is empty for distance from {pheno_from} to {pheno_to}")

    return df


#*****************************************************************

def calculate_median_distribution(dictionary_group,groups):
    dict_median={}
    for g in groups:
        dict_median[g]="NaN"
    
    for _k,_v in dictionary_group.items():
        if len (_v)==0:
            continue
        # supporta sia liste di float che liste di tuple (paz, valore)
        if isinstance(_v[0], tuple):
            only_values = [val for paz, val in _v]
        else:
            only_values = _v
        mean=np.median(only_values)
        dict_median[_k]=mean
        
    grade_major=""

    if all(value=="NaN"for value in dict_median.values()):
        grade_major="NaN"
        return grade_major,dict_median
    else:
        valid_value={k:v for k,v in dict_median.items() if str(v) !="NaN"}
        grade_major = max(valid_value, key=valid_value.get)
        return grade_major,dict_median




#*****************************************************************

def rearrange_lists(list_item):
    
    diff= abs(len(list_item[0])-len(list_item[1]))
    if diff!=0:
        list_len = [len(i) for i in list_item]
        i=[index for index, value in enumerate(list_len) if value == min(list_len)]
        list_item[i[0]].extend(np.full(shape=diff, fill_value=np.nan))
        
    return list_item
        
def statistical_test(df,path_output_result, test, p_adj):

    groups=list(df["GROUP"].unique())

    #array containing distance values for each cell and each subject
    values_distance=[]

    for g in groups:
        temp=list(df[df["GROUP"]==g]["DISTANCE"].values)
        values_distance.append(temp)

    p_value=10
    kruskal=False

    try:
        diff= abs(len(values_distance[0])-len(values_distance[1]))
    except Exception:
        diff= abs(len(values_distance[0]))
        
    if len(df["GROUP"].unique())==2 and diff !=0:
        logger.warning("The paired couple has a different observation number: Mann-Whitney test will be use!")
        test="unpaired"

    if len(df["GROUP"].unique())==1:
        logger.warning("Only One Group - not statistical is possible")
    
    #Case 2 groups---> Mann-Whitney test
    elif len(df["GROUP"].unique())==2 and test != "paired":
        logger.info("Running Mann-Whitney test")
        _, p_value = stats.mannwhitneyu(values_distance[0], values_distance[1], nan_policy='omit')
    
    #Case 2 groups---> Wilcoxon test (paired test)
    elif len(df["GROUP"].unique())==2 and test == "paired":
        #check list dimensions
        #values_distance=rearrange_lists(values_distance)
        logger.info("Running Wilcoxon test")
        _, p_value = stats.wilcoxon(values_distance[0], values_distance[1], nan_policy='omit')
    
    #Case more than 2 groups---> Kruskal test
    elif len(df["GROUP"].unique())>2:
        logger.info("Running Kruskal test")
        _, p_value = stats.kruskal(*[v for v in values_distance], nan_policy='omit')

    # Proceed with Dunn's test only if Kruskal-Wallis is significant   
        if p_value <= 0.05:  
            print(p_value)
            kruskal=True
            print("Running Dunn's test")
            if p_adj==None:
                p_adj="bonferroni"
            dunn_results = sp.posthoc_dunn([*[v for v in values_distance]],p_adjust=p_adj)
           # Convert results to DataFrame
            dunn_comparison_df = pd.DataFrame(dunn_results, columns=[i+1 for i in range(len(groups))], index=[i+1 for i in range(len(groups))])
            new_columns = {i: groups[i-1] for i in dunn_comparison_df.columns}
            new_index = {i: groups[i-1] for i in dunn_comparison_df.index}
            dunn_comparison_df=dunn_comparison_df.rename(columns=new_columns, index=new_index)
            dire =os.path.join(path_output_result,"distance_Dunn_Test")
            if not os.path.exists(dire):
                os.makedirs(dire)
            #dunn_comparison_df.to_csv(f"{dire}/Dunn_test_results.csv",sep="\t")
            dunn_comparison_df.to_csv(os.path.join(dire,"Dunn_test_results.csv"),sep="\t")
        
    return p_value,kruskal,test

#*****************************************************************

def rearrange_df(df):

    diff= abs(len(df[df["GROUP"]==df["GROUP"].unique()[0]])-len(df[df["GROUP"]==df["GROUP"].unique()[1]]))
    less_freq=df.GROUP.mode()[0]
    
    if diff!=0:
        df_new=pd.DataFrame(index=range(diff), columns=["GROUP", "DISTANCE"])
        df_new["GROUP"]=df[df["GROUP"]!=less_freq]["GROUP"].unique()[0]
        df=pd.concat([df, df_new])
        
    return df
def box_plots_distances(path_output_results, df, pheno_from, pheno_to, kruskal, p_adjust, test):

    pairwise_pvalues = {}  # <<<<<<<<<<<<<<<<<<<<<< AGGIUNTO

    if len(df["GROUP"].unique()) != 1:
        dire = os.path.join(path_output_results, "box_plot")
        os.makedirs(dire, exist_ok=True)
        filename = os.path.join(dire, f"distances_box_plot_{pheno_from}_to_{pheno_to}.png")

        groups = list(df["GROUP"].unique())

        # === 2 GRUPPI: MANN-WHITNEY O WILCOXON ===
        if len(groups) == 2:
            g1, g2 = groups
            x1 = df[df["GROUP"] == g1]["DISTANCE"].dropna().values
            x2 = df[df["GROUP"] == g2]["DISTANCE"].dropna().values

            if test != "paired":
                # Mann-Whitney esattamente come fa tap
                _, pval = stats.mannwhitneyu(x1, x2, nan_policy="omit")
            else:
                # Wilcoxon
                _, pval = stats.wilcoxon(x1, x2, nan_policy="omit")

            pairwise_pvalues[(g1, g2)] = pval

            # Plot originale
            tap.plot_stats(df,
                        x="GROUP",
                        y="DISTANCE",
                        type_correction=p_adjust,
                        filename=filename,
                        kwargs={"title":f'Distance from {pheno_from} to {pheno_to}',
                                "labels":{"GROUP": "Group",
                                            "DISTANCE": r'$Distance_{z}$',
                                            "GROUP": 'Group'}})



        # === 3+ GRUPPI: DUNN (stessi valori del boxplot) ===
        elif len(groups) > 2:

            if kruskal:
                if p_adjust is None:
                    p_adjust = "bonferroni"

                # Calcolo Dunn IDENTICO a tap
                value_lists = [df[df["GROUP"] == g]["DISTANCE"].dropna().values for g in groups]
                dunn_df = sp.posthoc_dunn(value_lists, p_adjust=p_adjust)

                # Rinomina righe/colonne con i nomi corretti dei gruppi
                dunn_df.index = groups
                dunn_df.columns = groups

                # Estrai i p-value ESATTI che vengono plottati
                for g1 in groups:
                    for g2 in groups:
                        if g1 != g2:
                            pairwise_pvalues[(g1, g2)] = dunn_df.loc[g1, g2]

                # Plot originale
                tap.plot_stats(
                    df, x="GROUP", y="DISTANCE",
                    type_test="dunn",
                    type_correction=p_adjust,
                    filename=filename,
                    kwargs={"title": f"Distance from {pheno_from} to {pheno_to}",
                            "labels": {"GROUP": "Group", "DISTANCE": r"$Distance_z$"}}
                )

            else:
                # Non significativo → solo il plot
                fig = px.box(df, x="GROUP", y="DISTANCE", color="GROUP",
                             title=f"Distance from {pheno_from} to {pheno_to} (Kruskal-Wallis p>0.05)")
                fig.update_traces(quartilemethod="linear")
                fig.write_image(filename)

    return pairwise_pvalues  # <<<<<<<<<<<<<<<<<<<<<< AGGIUNTO

#*****************************************************************

def plot_distance_curve(path_output_result,df,pheno_from,pheno_to,p_value):
    dire =os.path.join(path_output_result,"distance_curve")
    if not os.path.exists(dire):
        os.makedirs(dire)

    #group list
    groups=list(df["GROUP"].unique())

    #array with distance values for each cell and each subject
    values_distance=[]

    for g in groups:
        temp=list(df[df["GROUP"]==g]["DISTANCE"].values)
        values_distance.append(temp)
    
    fig = ff.create_distplot(values_distance, groups,show_hist=False,show_rug=False)
    filename=os.path.join(dire,"plot_statistical_distance_"+pheno_from+"_to_"+pheno_to+".png")

    if p_value < 0.05 and p_value >= 0.001:
        fig.update_layout({'plot_bgcolor':'white'},title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:green;'>p value < 0.05</span>",
                        yaxis=dict(tickformat=".4f",title_text=r"$Density$"),xaxis=dict(title_text=r"$Distance_{z}$"),legend_title_text="Group")
        fig.write_image(filename,scale=6)
    
    elif p_value < 0.001:
        fig.update_layout({'plot_bgcolor':'white'},title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:green;'>p value < 0.001</span>",
                        yaxis=dict(tickformat=".4f",title_text=r"$Density$",gridcolor='lightgrey'),xaxis=dict(title_text=r"$Distance_{z}$",gridcolor='lightgrey'),legend_title_text="Group")
        fig.write_image(filename,scale=6)
        
    elif p_value >= 0.05 and p_value<10:
        fig.update_layout({'plot_bgcolor':'white'},title_text=f"Distance-Z score from {pheno_from} to {pheno_to}<br> <span style ='font-size: 10px;color:red;'>p value > 0.05</span>",
                        yaxis=dict(tickformat=".4f",title_text=r"$Density$",gridcolor='lightgrey'),xaxis=dict(title_text=r"$Distance_{z}$",gridcolor='lightgrey'),legend_title_text="Group")
        fig.write_image(filename,scale=6)

#*****************************************************************
#*****************************************************************
def main(data):
    print("\n######################## STATISTICAL ANALYSIS DISTANCE #########################\n")
    logger.info("Start statistical distance analysis process")

    # path folder with distance
    root_folder = os.path.join(data["Paths"]["output_folder"], "Distance")
    if not os.listdir(root_folder):
        logger.critical(f"{root_folder} is an empty directory")
        return

    # path folder to save distance statistical results
    path_output = os.path.join(data["Paths"]["output_folder"], "Distance_Statistical")
    create_output_dir(path_output)
    logger.info(f"Statistical distances output will be stored in {path_output}")

    # list of phenotypes
    pheno_interested = data["Phenotypes"]["pheno_list"]
    pheno_from = data["Distance"]["pheno_from"]
    pheno_to = data["Distance"]["pheno_to"]
    plot_distance = data["Distance"]["plot_distance"]
    save_csv_zetascore = data["Distance"]["save_csv_zetascore"]

    p_adjust = data["Stats"]["p_adj"]
    if p_adjust is None or p_adjust.strip() == "":
        p_adjust = "bonferroni"
    else:
        p_adjust = p_adjust.lower()



    st_test = data["Stats"]["sample_type"]

    # groups
    groups = [f for f in os.listdir(root_folder) if not f.startswith('.')]
    logger.info(f"{len(groups)} group(s) found!")

    path_stats_file = os.path.join(path_output, "summary_statistical.csv")
    if os.path.exists(path_stats_file):
        os.remove(path_stats_file)

    df_all = pd.DataFrame(columns=["GROUP", "PHENO", "DISTANCE"])
    pvalues_dict = {}
    dict_statistical_result = {}

    # caso: tutti i fenotipi
    if pheno_from == "" and pheno_to == "":
        logger.info(f"No specific phenotype selected. Evaluating all combinations: {pheno_interested}")
        phenos_to_run = list(itertools.permutations(pheno_interested, 2))

        with ProcessPoolExecutor() as executor:
            futures = [
                executor.submit(process_pheno_pair, pheno, data, root_folder, path_output, groups,
                                plot_distance, st_test, p_adjust, save_csv_zetascore)
                for pheno in phenos_to_run
            ]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing pairs"):
                pheno_from_res, pheno_to_res, pvalue, grade_major, dict_median, pvalues = future.result()
                
                # Mantieni None se pvalue non calcolabile
                dict_statistical_result[f"{pheno_from_res}to{pheno_to_res}"] = {
                    "p_value": pvalue,  # None se non disponibile
                    "grade_major": grade_major,
                    "median": dict_median,
                }
                pvalues_dict.update(pvalues)

                df_temp = create_df_distances(
                    prepare_dataframe_distances(root_folder, pheno_from_res, pheno_to_res),
                    path_output, pheno_from_res, pheno_to_res, save_csv_zetascore
                )
                df_temp["PHENO"] = pheno_from_res
                df_temp["PHENO_FROM"] = pheno_from_res
                df_temp["PHENO_TO"] = pheno_to_res
                df_all = pd.concat([df_all, df_temp], ignore_index=True)

    else:
        logger.info(f"pheno_from={pheno_from}, pheno_to={pheno_to} selected")
        if pheno_from not in pheno_interested or pheno_to not in pheno_interested:
            logger.critical("Invalid pheno_from or pheno_to")
            return

        pheno_result = process_pheno_pair((pheno_from, pheno_to), data, root_folder, path_output,
                                          groups, plot_distance, st_test, p_adjust, save_csv_zetascore)
        pheno_from_res, pheno_to_res, pvalue, grade_major, dict_median, pvalues = pheno_result

        dict_statistical_result[f"{pheno_from_res}to{pheno_to_res}"] = {
            "p_value": pvalue,  # None se non calcolabile
            "grade_major": grade_major,
            "median": dict_median,
        }
        pvalues_dict.update(pvalues)

        df_all = create_df_distances(
            prepare_dataframe_distances(root_folder, pheno_from_res, pheno_to_res),
            path_output, pheno_from_res, pheno_to_res, save_csv_zetascore
        )
        df_all["PHENO"] = pheno_from_res

    # salva file riassuntivo
    # salva file riassuntivo
    with open(path_stats_file, "w") as f:
        f.write("Phenotype\tDistance_to\tP_value\t" + "\t".join(f"Median_{g}" for g in groups) + "\tGroup_major\n")
        for pheno, stat in dict_statistical_result.items():
            ph = pheno.split("to")
            # P-value
            pval_str = "None" if stat["p_value"] is None else stat["p_value"]
            # Mediane
            median_strs = []
            median_dict = stat.get("median") or {}
            for g in groups:
                val = median_dict.get(g)
                if val is None or (isinstance(val, str) and val.lower() == "nan"):
                    median_strs.append("None")
                else:
                    median_strs.append(val)

            f.write(f"{ph[0]}\t{ph[1]}\t{pval_str}\t" + "\t".join(map(str, median_strs)) + f"\t{stat['grade_major']}\n")


    # heatmap NxN
    if not df_all.empty:
        plot_heatmap_and_bar(df_all, path_output, phenotypes=pheno_interested, pvalues_dict=pvalues_dict)

    merge_all_pairwise_csv(path_output)
    logger.info("End statistical analysis!")
