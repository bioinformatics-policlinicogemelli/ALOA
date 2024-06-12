import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import os
import numpy as np
from image_proc_functions import pheno_filt
import json
import pathlib
import random
from datetime import datetime
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import StandardScaler
from loguru import logger

from kmodes.kprototypes import KPrototypes
from kneed import KneeLocator

import warnings
warnings.filterwarnings("ignore")

def plot_elbow_analysis(df, k_number, idx, output_path):
    logger.info("Creating the output folder " + os.path.join(output_path, "elbow_scores"))
    os.makedirs(os.path.join(output_path, "elbow_scores"), exist_ok=True)

    df_sample = df[['Cell.X.Position' , 'Cell.Y.Position']]

    logger.info("Starting elbow analysis")
    inertia_values = []
    K_range = range(2, k_number)

    for k in K_range:
        #logger.debug(k)
        print(f"Calcolo del grafico per k={k}")
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(df_sample.values)
        inertia_values.append(kmeans.inertia_)

    df_sample["Pheno"]=df["Pheno"]
    knee_locator = KneeLocator(list(K_range), inertia_values, curve="convex", direction="decreasing")
    n_opt_cluster = knee_locator.knee
    knee_locator.plot_knee()

    logger.info("Creating optimal cluster number graph") 
    # Visualizing inertia plot
    plt.plot(K_range, inertia_values, marker='o')
    plt.xlabel('N Cluster')
    plt.ylabel('Inertia')
    plt.title('Elbow method to find optimal number of cluster')
    plt.savefig(os.path.join(os.path.join(output_path, "elbow_scores"),idx+".tiff"), dpi=300, bbox_inches = "tight")
    plt.close()

    logger.info(f"Optimal cluster Number: {n_opt_cluster}")
    
    return df_sample, n_opt_cluster

def plot_silhouette_analysis(df, output_path, idx, k_number):

    logger.info("Creating the output folder " + os.path.join(output_path, "silhoutte_scores"))
    os.makedirs(os.path.join(output_path, "silhoutte_scores"), exist_ok=True)

    df_sample = df[['Cell.X.Position', 'Cell.Y.Position']]

    logger.info("Starting silhouette analysis")
    silhouette_scores = []
    etichette_cluster={}
    K_range = range(2, k_number)
    for k in K_range:
        logger.info(f"Check for {k} clusters")
        spectral = SpectralClustering(n_clusters=k, affinity='nearest_neighbors', random_state=42)
        #kmeans = KMeans(n_clusters=k)
        etichette_cluster[k] = spectral.fit_predict(df_sample.values)
        silhouette_avg = silhouette_score(df_sample.values, etichette_cluster[k])
        silhouette_scores.append(silhouette_avg)

    df_sample["Pheno"]=df["Pheno"]
    logger.info("Creating optimal cluster number graph")
    plt.plot(K_range, silhouette_scores, marker='o')
    plt.xlabel('N cluster')
    plt.ylabel('Silhouette score')
    plt.savefig(os.path.join(os.path.join(output_path, "silhoutte_scores"),idx+".tiff"), dpi=300, bbox_inches = "tight")
    plt.close()

    # Select optimal cluster number based on silhouette score
    n_opt_cluster = K_range[silhouette_scores.index(max(silhouette_scores))]
    logger.info(f"Optimal cluster Number: {n_opt_cluster}")
    
    return df_sample, n_opt_cluster

def clustering_function (n_opt_cluster, df_sample, clust_alg, df):
    
    if not len(clust_alg.lower().replace("k","").replace("s","").replace("p",""))==0:
        logger.critical(f"No valid clustering algorithm selected. Clust_alg '{clust_alg}' were selected but only 'kmeans' (k), 'spectral' (s) and 'prototype' (p) are currently supported! Check your cluster section in your config.json file!")
        exit()
        
    if "k" in clust_alg.lower():
        km = KMeans(n_clusters=n_opt_cluster)
        logger.info("kmeans algorithm selected!")
        df_sample['kmeans'] = km.fit_predict(df_sample[["Cell.X.Position", "Cell.Y.Position"]].values)

    if "s" in clust_alg.lower():
        logger.info("spectral algorithm selected!")
        
        scaler= StandardScaler()
        spatial_data_normalizer = scaler.fit_transform(df_sample[["Cell.X.Position", "Cell.Y.Position"]].values)
        spectral = SpectralClustering(n_clusters=n_opt_cluster, affinity='nearest_neighbors', random_state=42)
        df_sample['spectral']= spectral.fit_predict(spatial_data_normalizer)
    
    if "p" in clust_alg.lower():
        logger.info("k-prototype algoritm selected!")
        cat_cols= [2]
        kproto = KPrototypes(n_clusters=n_opt_cluster, init='Cao', verbose=2)
        df_sample['prototype'] = kproto.fit_predict(df_sample, categorical=cat_cols)

    df_sample["Pheno"]= df["Pheno"]
    
    return n_opt_cluster, df_sample

def plot_convex_hull(df_plot, n_opt_cluster, idx, output_path, clust_alg):
   
    output_path=os.path.join(output_path,clust_alg)
    os.makedirs(output_path, exist_ok=True)
        
    output_f=os.path.join(output_path, "scatter_plot")
    os.makedirs(output_f, exist_ok=True)

    phenolist=list(df_plot["Pheno"].unique())

    _, _ = plt.subplots()

    # color map for cluster
    cmap = plt.get_cmap("tab10")
    colors = cmap(np.arange(n_opt_cluster))
 
    logger.info("Creation of convex hull plot")
    for cluster in range(n_opt_cluster):

        cluster_data = df_plot[df_plot[clust_alg] == cluster]
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
        logger.error("Error: Empty DataFrame")
        return
    
    # Check for cluster column existance
    if clust_alg not in df_plot.columns:
        logger.error("Error: cluster method column not found in DataFrame")
        return

    # Check for pheno column existence
    if 'Pheno' not in df_plot.columns:
        logger.error("Error: 'Pheno' column not found in DataFrame")
        return
    
    output_barplot=os.path.join(output_path, "stacked_barplot")
    os.makedirs(output_barplot, exist_ok=True)

    df_stacked = df_plot.groupby([clust_alg, 'Pheno']).size().unstack(fill_value=0)

    # Percentage evaluation for each phenotype in each cluster
    cluster_percentages = df_stacked.div(df_stacked.sum(axis=1), axis=0)

    # stacked barplot creation
    fig,ax=plt.subplots()

    for i in range(df_stacked.values.shape[1]):
        ax.bar(
            sorted(list(df_plot[clust_alg].unique())),    
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
    
    print("\n############################### CLUSTER ANALYSIS ###############################\n")
    
    logger.info("Start clustering process: This step will produce a spatial clustering\n")

    logger.info("Reading configuration file")
    with open("config.json") as f:
        data=json.load(f)
        
    input_path = os.path.join(data["Paths"]["output_folder"],"Merged_clean")
   
    output_path = os.path.join(data["Paths"]["output_folder"],"Clustering")
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
    logger.info(f"Clustering output will be stored in {output_path}")
    
    clust_alg=data["cluster"]["cluster_method"]

    if not data["cluster"]["pheno_list"]:
        pheno_list=data["Phenotypes"]["pheno_list"]
    else:
        pheno_list=data["cluster"]["pheno_list"]
    
    groups=[f for f in os.listdir(input_path) if not f.startswith('.')] 

    for g in groups:

        output_path_g=os.path.join(output_path,g)
        pzt_list=[f for f in os.listdir(os.path.join(input_path,g)) if not f.startswith('.')] 
        logger.info("Starting analyzing group: "+g)
        
        for pzt in pzt_list:
            
            print("sbj: "+pzt)
            
            # Load and filter data
            df = pd.read_csv(os.path.join(input_path,g,pzt), sep='\t')
            df.columns = df.columns.str.replace(' ', '.')
            df_filt = pheno_filt(df, pheno_list)
            if len(df_filt)==0:
                print("No Phenotype(s) found. Check the phenotype list and Phenotype columns of your data!")
                print("Skip to next subject")
                continue

            idx=df["Sample.Name"][0].split(" ")[0]
            
            logger.warning(f"The ideal pairing for silhouette analysis is spectral clustering, for prototype analysis it's elbow method, and k-means is suitable for both.")
            
            number_alg=data["cluster"]["algo_method"].lower()
            k=data["cluster"]["k"]
            
            if number_alg=="e":
                #elbow method
                df_sample, n_cluster_opt= plot_elbow_analysis(df_filt, k, idx, output_path_g)
            if number_alg =="s":
                # silhouette analysis
                df_sample, n_cluster_opt = plot_silhouette_analysis(df_filt, output_path_g, idx, k)

            if number_alg not in ('e', 's'):
                    logger.critical(f"ERROR: wrong clustering algoritm for optimal number of cluster. '{number_alg}' were selected but only 'elbow method' (e) or 'silhouette analysis' (s) are currently supported! Check your method section in your config.json file!")
                    exit(1)

            _, df_sample= clustering_function (n_cluster_opt, df_sample, clust_alg, df)

            cluster_type=list(clust_alg)
            
            for cl in cluster_type:

                    if cl=="k":
                        cl_name="kmeans"
                    if cl=="s":
                        cl_name="spectral"
                    if cl=="p":
                        cl_name="prototype"
                    
                    # Convex Hull plot
                    output_res=plot_convex_hull(df_sample, n_cluster_opt, idx, output_path_g, cl_name)

                    # stacked barplot and percentage csv 
                    plot_stacked_bar_chart(df_sample, output_res, idx, cl_name)
    
    logger.info("End clustering analysis!\n")
    return ()

if __name__ == "__main__":
    main()