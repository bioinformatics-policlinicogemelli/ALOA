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
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import StandardScaler

import warnings
warnings.filterwarnings("ignore")

def plot_silhouette_analysis(df, output_path, idx, k_number, clust_alg):

    print("Creating the output folder " + os.path.join(output_path, "silhoutte_scores"))
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
    
    if clust_alg.lower()=="kmeans":
        df_sample['Cluster'] = km.fit_predict(df_sample.values)

    elif clust_alg.lower()=="spectral":
        scaler= StandardScaler()
        spatial_data_normalizer = scaler.fit_transform(df_sample[["Cell.X.Position", "Cell.Y.Position"]].values)
        spectral = SpectralClustering(n_clusters=n_opt_cluster, affinity='nearest_neighbors', random_state=42)
        df_sample['Cluster']= spectral.fit_predict(spatial_data_normalizer)
    
    else:
        print(f"ERROR: wrong clustering algorithm selected. '{clust_alg}' were selected but only 'kmeans' and 'spectral' are currently supported! Check your cluster section in your config.json file!")
        exit(1)

    df_sample["Pheno"]= df["Pheno"]
    
    return df_sample, n_opt_cluster, clust_alg

def plot_convex_hull(df_plot, n_opt_cluster, idx, output_path, clust_alg):

    output_path=os.path.join(output_path,clust_alg)
    os.makedirs(output_path, exist_ok=True)
        
    output_f=os.path.join(output_path, "scatter_plot")
    os.makedirs(output_f, exist_ok=True)

    phenolist=list(df_plot["Pheno"].unique())

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

        markers=['.','o','v','^','<','>',
                 '1','2','3','4',
                 's','p','*','h',
                 '+','x','D','d']
        random.seed(42)
        random.shuffle(markers)

        # plot for each phenotype
        for i,ph in enumerate(phenolist):
            points_ch = cluster_data[cluster_data['Pheno'].str.strip() == ph]
            if points_ch.empty:
                continue
            plt.scatter(points_ch['Cell.X.Position'], points_ch['Cell.Y.Position'], marker=markers[i], color=colors[cluster], label=f'C{cluster + 1} {ph}', s=4)
        
        for simplex in hull.simplices:
            plt.fill(points[simplex, 0], points[simplex, 1], edgecolor=colors[cluster], facecolor=colors[cluster], alpha=1, closed=True)

    plt.xlabel('Cell X Position')
    plt.ylabel('Cell Y Position')
    plt.title('Convex Hull')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(os.path.join(output_f,idx+".tiff"), dpi=300, bbox_inches = "tight")
    plt.close()
    return output_path

def plot_stacked_bar_chart(df_plot, output_path, idx, clust_alg):
    
    # Check for empty dataframe
    if df_plot.empty:
        print("Error: Empty DataFrame")
        return
    
    # Check for cluster column existance
    if 'Cluster' not in df_plot.columns:
        print("Error: 'Cluster' column not found in DataFrame")
        return

    # Check for pheno column existence
    if 'Pheno' not in df_plot.columns:
        print("Error: 'Pheno' column not found in DataFrame")
        return
    
    output_barplot=os.path.join(output_path, "stacked_barplot")
    os.makedirs(output_barplot, exist_ok=True)

    df_stacked = df_plot.groupby(['Cluster', 'Pheno']).size().unstack(fill_value=0)

    # Calcolo delle percentuali per ogni Pheno nei vari cluster
    cluster_percentages = df_stacked.div(df_stacked.sum(axis=1), axis=0)

    # stacked barplot creation
    fig,ax=plt.subplots()

    for i in range(df_stacked.values.shape[1]):
        ax.bar(
            sorted(list(df_plot["Cluster"].unique())),    
            df_stacked.values.transpose()[i], bottom = np.sum(df_stacked.values.transpose()[:i], axis = 0), width = 0.25
            )
    
    ax.set_xticks(list(range(0, df_stacked.values.shape[0])))
    ax.set_xticklabels(list(range(0, df_stacked.values.shape[0])))
    ax.legend(labels=list(df_stacked.columns), title='Phenotype', loc='upper left', bbox_to_anchor=(1, 1), frameon=False)

    plt.xlabel('Cluster')
    plt.ylabel('Count')
    plt.title('Phenotype(s) distribution inside cluster')
    plt.tight_layout() 
    plt.savefig(os.path.join(output_barplot,idx+".tiff"), dpi=300, bbox_inches = "tight")
    plt.close()

    output_barplot=os.path.join(output_path, "percentage")
    os.makedirs(output_barplot, exist_ok=True)

    # plot phenotype percentage inside cluster
    cluster_percentages=cluster_percentages * 100
    cluster_percentages.to_csv(os.path.join(output_barplot,"cluster_percentage_"+idx+".csv"), sep=",")

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
            clust_alg=data["cluster"]["cluster_method"]
            df_sample_tmp, n_cluster_opt, clust_alg = plot_silhouette_analysis(df_filt, output_path_g, idx, k, clust_alg)

            # Convex Hull plot
            output_res=plot_convex_hull(df_sample_tmp, n_cluster_opt, idx, output_path_g, clust_alg)

            # stacked barplot and percentage csv 
            plot_stacked_bar_chart(df_sample_tmp, output_res, idx, clust_alg)

if __name__ == "__main__":
    main()
