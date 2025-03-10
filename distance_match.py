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

from image_proc_functions import pheno_filt, img_filt, crop_img, norm_values, plot_pheno, plot_dist, dist_eval
import cv2
import pandas as pd
import os
import pathlib
import itertools
import re
from loguru import logger
from pandas.api.types import is_numeric_dtype

def distance_match(data):

    print("\n################################ DISTANCE MATCH ################################\n")
    
    logger.info("Start distance match process: This step will provide the plot of the distance between couple of phenotypes on the tiff image(s)\n")
    
    output=data["Paths"]["output_folder"]
    output_f=os.path.join(output,"Distance_match")
    pathlib.Path(output_f).mkdir(parents=True, exist_ok=True)

    img_folder=os.path.join(data["Paths"]["data_input_folder"],"img_match")
    if not os.path.exists(img_folder):
        logger.error(f"No image folder found in input folder!") 
        return 
    csv_folder=os.path.join(data["Paths"]["data_input_folder"],"raw_data") 
    
    if len(data["Distance_match"]["pheno_list"])==0:
        pheno_list=data["Phenotypes"]["pheno_list"]
    else:
        pheno_list=data["Distance_match"]["pheno_list"]
    
    dir=[f for f in os.listdir(img_folder) if not f.startswith('.')]
    img_path=list(set(list(map(lambda x: os.path.join(img_folder,x), dir))))
    csv_path=list(set(list(map(lambda x: os.path.join(csv_folder,x.split("_[")[0],re.search('.*?\]', x).group(0) +"_cell_seg_data.txt"), dir))))
    
    img_path.sort()
    csv_path.sort()
    
    for i,c in zip(img_path,csv_path):
        
        print(i)
        print(c)

        logger.info(f"Analyzing file {c}")
        
        #dataframe section
        try:
            df=pd.read_csv(c,sep="\t")
            df.columns=list(map(lambda x: x.replace(" ",".") ,df.columns))
            if not is_numeric_dtype(df["Cell.X.Position"]):
                df=pd.read_csv(c, sep="\t", decimal=",")
                df.columns=list(map(lambda x: x.replace(" ",".") ,df.columns))
        except FileNotFoundError:
            logger.warning(f"No image matched file found in {c}. Skip to the next image!")
            continue
        
        df.columns = df.columns.str.replace(' ', '.')
        
        # image section
        img = cv2.imread(i)
        masked_img=img_filt(img)
        crop=crop_img(masked_img,50)

        if crop.shape[0]==img.shape[0] and crop.shape[1]==img.shape[1]: 
            crop=img
        else:
            if crop.shape[0]<img.shape[0]:
                pass
            if crop.shape[1]<img.shape[1]:
                pass
        
        #Normalize Data
        pheno_df=norm_values(df, crop)
        comb=list(itertools.combinations(pheno_list, 2))
        for cb in comb:
            logger.info(f"Check distance between {cb[0]} and {cb[1]}")
            
            pheno_df=pheno_filt(df,list(cb))
            if len(pheno_df["Pheno"].unique())!=2:
                logger.warning(f"Data {c} has no match with requested Phenotype list. Skip to the next one!")
                continue

            pheno_df_dist = dist_eval(pheno_df)
            
            #plot
            ph1=pheno_df_dist[pheno_df_dist["Phenotype"]==cb[0]]
            ph2=pheno_df_dist[pheno_df_dist["Phenotype"]==cb[1]]
            
            if len(ph1)==0:
                logger.warning(f"{ph1} not found!")
                logger.info("skip to the next couple oh phenotypes")
                continue
            
            if len(ph2)==0:
                logger.warning(f"{ph2} not found!")
                logger.info("skip to the next couple oh phenotypes")
                continue
            
            filename=os.path.splitext(os.path.basename(i))[0]+"_dist_match_" + ''.join(map(str, list(cb)))
            #filename=i.split("/")[-1].split("/")[-1].replace(os.path.splitext(i)[-1],"_dist_match_") + ''.join(map(str, list(cb)))
            
            # nearest ph2 cell for each ph1
            ph1_to_ph2=pd.merge(ph1, ph2, left_on='Cell ID ' + cb[1], right_on='Cell ID', suffixes=[None, "."+cb[1]])
            if len(ph1_to_ph2)>0:
                plot_dist(crop, pheno_df_dist, ph1_to_ph2, os.path.join(output_f, filename))
            
            # nearest ph1 cell for each ph2
            ph2_to_ph1=pd.merge(ph2, ph1, left_on='Cell ID ' + cb[0], right_on='Cell ID', suffixes=[None, "."+cb[0]])
            if len(ph2_to_ph1)>0:
                plot_dist(crop, pheno_df_dist, ph2_to_ph1, os.path.join(output_f, filename))
            
            # mutual nearest neighbors (cells which have each other as nearest neighbors)
            mutual=ph2_to_ph1[ph2_to_ph1["Cell ID"]==ph2_to_ph1["Cell ID "+cb[1]+"."+cb[0]]]
            if len(mutual)>0:
                plot_dist(crop, pheno_df_dist, mutual, os.path.join(output_f, filename),m=True)
    
    logger.info("End distance match!\n")
    return()

if __name__ == '__main__':
    distance_match()