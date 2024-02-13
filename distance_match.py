#libraries import
from image_proc_functions import pheno_filt, img_filt, crop_img, norm_values, plot_pheno, plot_dist, dist_eval
import cv2
import pandas as pd
import json
import os
import pathlib
import itertools

def distance_match():
    
    with open('config.json', 'r') as myfile:
        data=myfile.read()
    # parse json file
    obj = json.loads(data)

    output=obj["Paths"]["output_folder"]
    pathlib.Path(output).mkdir(parents=True, exist_ok=True) 
    output_f=os.path.join(output,"Distance_match")
    pathlib.Path(output_f).mkdir(parents=True, exist_ok=True)

    input_folder=os.path.join(obj["Paths"]["data_input_folder"],"img_match")
    
    if len(obj["distance_match"]["pheno_list"])==0:
        pheno_list=obj["Phenotypes"]["pheno_list"]
    else:
        pheno_list=obj["distance_match"]["pheno_list"]
    
    dir=os.listdir(input_folder)
    dir=list(set(list(map(lambda x: os.path.join(input_folder,x.split(']_')[0]+"]"), dir))))

    for d in dir:
        print(f"Analyzing file {d}")
        #dataframe section
        df=pd.read_csv(d+"_cell_seg_data.txt",sep="\t")
        
        # image preproc section
        img = cv2.imread(d+"_composite_image.tif")
        masked_img=img_filt(img)
        crop=crop_img(masked_img,50)

        if crop.shape[0]==img.shape[0] and crop.shape[1]==img.shape[1]: 
            crop=img
        else:
            if crop.shape[0]<img.shape[0]:
                pass
            if crop.shape[1]<img.shape[1]:
                pass
        
        pheno_df=norm_values(df, crop)
        comb=list(itertools.combinations(pheno_list, 2))
        for c in comb:
            print(f"Check distance between {c[0]} and {c[1]}")
            
            pheno_df=pheno_filt(df,list(c))
            if len(pheno_df["Pheno"].unique())!=2:
                print("Data " + d + "has no match with requested Phenotype list. Skip to the next one!")
                continue

            
            pheno_df_dist = dist_eval(pheno_df)
            
            #plot
            ph1=pheno_df_dist[pheno_df_dist["Phenotype"]==c[0]]
            ph2=pheno_df_dist[pheno_df_dist["Phenotype"]==c[1]]
            
            if len(ph1)==0:
                print(f"{ph1} not found!")
                print("skip to the next couple oh phenotypes")
                continue
            
            if len(ph2)==0:
                print(f"{ph2} not found!")
                print("skip to the next couple oh phenotypes")
                continue
            
            filename=(d+"_composite_image.tif").split("/")[-1].replace(".tif","_dist_match_") + ''.join(map(str, list(c)))
            
            # nearest ph2 cell for each ph1
            ph1_to_ph2=pd.merge(ph1, ph2, left_on='Cell ID ' + c[1], right_on='Cell ID', suffixes=[None, "."+c[1]])
            if len(ph1_to_ph2)>0:
                plot_dist(crop, pheno_df_dist, ph1_to_ph2, os.path.join(output_f, filename))
            
            # nearest ph1 cell for each ph2
            ph2_to_ph1=pd.merge(ph2, ph1, left_on='Cell ID ' + c[0], right_on='Cell ID', suffixes=[None, "."+c[0]])
            if len(ph2_to_ph1)>0:
                plot_dist(crop, pheno_df_dist, ph2_to_ph1, os.path.join(output_f, filename))
            
            # mutual nearest neighbors (cells which have each other as nearest neighbors)
            mutual=ph2_to_ph1[ph2_to_ph1["Cell ID"]==ph2_to_ph1["Cell ID "+c[1]+"."+c[0]]]
            if len(mutual)>0:
                plot_dist(crop, pheno_df_dist, mutual, os.path.join(output_f, filename),m=True)

if __name__ == '__main__':
    distance_match()