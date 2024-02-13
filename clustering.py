import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import os
import seaborn as sns
import numpy as np
from image_proc_functions import pheno_filt
import json
import pathlib
import random

import warnings
warnings.filterwarnings("ignore")

def plot_silhouette_analysis(df, output_path, idx, k_number):

    print("Creating the output folder "+os.path.join(output_path, "silhoutte_scores"))
    os.makedirs(os.path.join(output_path, "silhoutte_scores"), exist_ok=True)

    df_sample = df[['Cell.X.Position', 'Cell.Y.Position']]

    print("silhouette evaluation")
    silhouette_scores = []
    K_range = range(2, k_number)
    for k in K_range:
        print(k)
        kmeans = KMeans(n_clusters=k)
        etichette_cluster = kmeans.fit_predict(df_sample.values)
        silhouette_avg = silhouette_score(df_sample.values, etichette_cluster)
        silhouette_scores.append(silhouette_avg)

    print("Creating optimal cluster number graph")
    plt.plot(K_range, silhouette_scores, marker='o')
    plt.xlabel('N cluster')
    plt.ylabel('Silhouette score')
    plt.savefig(os.path.join(os.path.join(output_path, "silhoutte_scores"),idx+".tiff"), dpi=300, bbox_inches = "tight")
    plt.close()

    # Select optimal cluster number based on silhouette score
    n_opt_cluster = K_range[silhouette_scores.index(max(silhouette_scores))]
    print(f"Optimal cluster Number: {n_opt_cluster}")

    km = KMeans(n_clusters=n_opt_cluster)
    label = km.fit_predict(df_sample.values)
    print(df_sample.head())  # O stampa df_sample.columns

    df_sample["Pheno"]= df["Pheno"]

    return df_sample, label, n_opt_cluster

def plot_convex_hull(df_plot, label, n_opt_cluster, idx, output_path):

    os.makedirs(os.path.join(output_path, "scatter_plot"), exist_ok=True)

    phenolist=list(df_plot["Pheno"].unique())
    df_plot['Cluster'] = label

    fig, ax = plt.subplots()

    # color map for cluster
    cmap = plt.get_cmap("tab10")
    colors = cmap(np.arange(n_opt_cluster))
 
    print("Creation of convex hull plot")
    for cluster in range(n_opt_cluster):

        cluster_data = df_plot[df_plot['Cluster'] == cluster]

        points = cluster_data[['Cell.X.Position', 'Cell.Y.Position']].values

        # convex hull evaluation
        hull = ConvexHull(points)

        for simplex in hull.simplices:
            plt.fill(points[simplex, 0], points[simplex, 1], edgecolor=colors[cluster], facecolor=colors[cluster], alpha=1, closed=True)

        #list of markers
        markers=['.','o','v','^','<','>',
                 '1','2','3','4',
                 's','p','*','h',
                 '+','x','D','d']
        random.seed(42)
        random.shuffle(markers)

        # plot for each phenotype
        for i,ph in enumerate(phenolist):
            points = cluster_data[cluster_data['Pheno'].str.strip() == ph]
            if points.empty:
                continue
            plt.scatter(points['Cell.X.Position'], points['Cell.Y.Position'], marker=markers[i], color=colors[cluster], label=f'C{cluster + 1} {ph}', s=2.5)

    plt.xlabel('Cell X Position')
    plt.ylabel('Cell Y Position')
    plt.title('Risultati con Convex Hull Colorati')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(os.path.join(os.path.join(output_path, "scatter_plot"),idx+".tiff"), dpi=300, bbox_inches = "tight")
    plt.close()

def plot_stacked_bar_chart(df_plot, output_path, idx):
    # Verifica se la colonna 'Cluster' esiste nel DataFrame

    if 'Cluster' not in df_plot.columns:
        print("Errore: La colonna 'Cluster' non è presente nel DataFrame.")
        return

    # Verifica se la colonna 'Pheno' esiste nel DataFrame
    if 'Pheno' not in df_plot.columns:
        print("Errore: La colonna 'Pheno' non è presente nel DataFrame.")
        return

    # Verifica se il DataFrame non è vuoto
    if df_plot.empty:
        print("Errore: Il DataFrame è vuoto.")
        return
    
    os.makedirs(os.path.join(output_path, "stacked_barplot"), exist_ok=True)

    df_stacked = df_plot.groupby(['Cluster', 'Pheno']).size().unstack(fill_value=0)

    # Calcolo delle percentuali per ogni Pheno nel totale
    total_percentages = df_plot['Pheno'].value_counts(normalize=True)

    # Calcolo delle percentuali per ogni Pheno nei vari cluster
    cluster_percentages = df_stacked.div(df_stacked.sum(axis=1), axis=0)

    # Ottieni i colori associati ai valori unici nella colonna 'Pheno'
    # palette = sns.color_palette('Set1')
    # unique_pheno_colors = {pheno: palette[i % len(palette)] for i, pheno in enumerate(df_plot['Pheno'].unique())}
    # legend_labels = [plt.Line2D([0], [0], marker='o', color='w', label=f'{pheno}', markerfacecolor=unique_pheno_colors[pheno], markersize=10) for pheno in df_plot['Pheno'].unique()]

    # stacked barplot creation
    fig,ax=plt.subplots()

    for i in range(df_stacked.values.shape[1]):
        ax.bar(
            sorted(list(df_plot["Cluster"].unique())),    
            df_stacked.values.transpose()[i], bottom = np.sum(df_stacked.values.transpose()[:i], axis = 0), width = 0.25
            )
    
    #ax.legend(loc="upper right")
    ax.set_xticks(list(range(0, df_stacked.values.shape[0])))
    ax.set_xticklabels(list(range(0, df_stacked.values.shape[0])))
    ax.legend(labels=list(df_stacked.columns), title='Phenotype', loc='upper left', bbox_to_anchor=(1, 1), frameon=False)

    plt.xlabel('Cluster')
    plt.ylabel('Count')
    plt.title('Distribuzione di Pheno nei vari Cluster')
    plt.tight_layout() 
    plt.savefig(os.path.join(os.path.join(output_path, "stacked_barplot"),idx+".tiff"), dpi=300, bbox_inches = "tight")
    #plt.show()
    plt.close()

    os.makedirs(os.path.join(output_path, "percentage"), exist_ok=True)

    # save percentage csv
    total_percentages=total_percentages * 100
    total_percentages.to_csv(os.path.join(os.path.join(output_path, "percentage"),"total_percentage_"+idx+".csv"), sep=",")

    # Stampa delle percentuali nei vari cluster
    cluster_percentages=cluster_percentages * 100
    cluster_percentages.to_csv(os.path.join(os.path.join(output_path, "percentage"),"cluster_percentage_"+idx+".csv"), sep=",")

def main():
    
    f=open("config.json")
    data=json.load(f)
    
    if data["Clean_data"]["other_rm"]:
        input_path = os.path.join(data["Paths"]["output_folder"],"Merged_clean")
    else:
        input_path = os.path.join(data["Paths"]["output_folder"],"Merged")
        
    output_path = os.path.join(data["Paths"]["output_folder"],"Clustering")
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)  
    
    groups=[f for f in os.listdir(input_path) if not f.startswith('.')] 

    if not data["cluster"]["pheno_list"]:
        pheno_list=data["Phenotypes"]["pheno_list"]
    else:
        pheno_list=data["cluster"]["pheno_list"]
    
    for g in groups:

        output_path_g=os.path.join(output_path,g)
        pzt_list=[f for f in os.listdir(os.path.join(input_path,g)) if not f.startswith('.')] 
        print("Starting analyzing group: "+g)
        for pzt in pzt_list:
            
            print("sbj: "+pzt)
            # Load and filter data
            df = pd.read_csv(os.path.join(input_path,g,pzt), sep='\t')
            df_filt = pheno_filt(df, pheno_list)
            idx=df["Sample.Name"][0].split(" ")[0]
            
            # silhouette analysis
            k=data["cluster"]["k"]
            df_sample_tmp, cluster_label, n_cluster_opt = plot_silhouette_analysis(df_filt, output_path_g, idx, k)

            # Convex Hull plot
            plot_convex_hull(df_sample_tmp, cluster_label, n_cluster_opt, idx, output_path_g)

            # stacked barplot and percentage csv 
            plot_stacked_bar_chart(df_sample_tmp, output_path_g, idx)

if __name__ == "__main__":
    main()
