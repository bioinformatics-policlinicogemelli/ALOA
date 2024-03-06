#libraries import
from image_proc_functions import pheno_filt, img_filt, crop_img, norm_values, plot_pheno, plot_interactive
import cv2
import pandas as pd
import json
import os
import pathlib
from loguru import logger

def img_match():
    
    logger.info("############################# DESCRIPTIVE ANALYSIS #############################")
    
    with open("config.json") as f:
        data=json.load(f)
    
    output=data["Paths"]["output_folder"]
            
    input_folder=os.path.join(data["Paths"]["data_input_folder"],"img_match") 
    if len(data["distance_match"]["pheno_list"])==0:
        pheno_list=data["Phenotypes"]["pheno_list"]
    else:
        pheno_list=data["image_match"]["pheno_list"]

    output_f=os.path.join(output,"Img_match")
    pathlib.Path(output_f).mkdir(parents=True, exist_ok=True)
    
    dir=[f for f in os.listdir(input_folder) if not f.startswith('.')]
    dir=list(set(list(map(lambda x: os.path.join(input_folder,x.split(']_')[0]+"]"), dir))))
    
    for d in dir:
        
        print("Analyzing data "+ d)
        
        #dataframe section
        df=pd.read_csv(d+"_cell_seg_data.txt",sep="\t")
        
        # image 
        img_filename=d+"_composite_image"
        if os.path.isfile(img_filename+".tif"):
            img_filename=img_filename+".tif"
        elif os.path.isfile(img_filename+".png"):
            img_filename=img_filename+".tif"
        elif os.path.isfile(img_filename+".jpg"):
            img_filename=img_filename+".jpg"
        else:
            logger.critical("Image not found! Exiting img match script!")
        img = cv2.imread(img_filename)
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
        
        #Filter by phenotype of interest
        pheno_df=pheno_filt(pheno_df,pheno_list)
        
        if len(pheno_df)==0:
            print("Data " + d + "has no match with requested Phenotype list. Skip to the next one!")
            continue
        
        #plot
        filename=(d+"_composite_image.tif").split("/")[-1].replace(".tif","_match_") + ''.join(map(str, pheno_list))
        plot_pheno(crop, pheno_df, os.path.join(output_f, filename))
        
        if data["image_match"]["interactive"]:
            pathlib.Path(os.path.join(output_f, "Interactive_plots")).mkdir(parents=True, exist_ok=True) 
            plot_interactive(crop, pheno_df, os.path.join(os.path.join(output_f, "Interactive_plots"), filename), 
                             data["image_match"]["layout_marker_edge_col"], 
                             data["image_match"]["layout_marker_size"], 
                             data["image_match"]["layout_xsize"], data["image_match"]["layout_ysize"])
        
    print("All done!")
    
if __name__ == '__main__':
    img_match()