from image_proc_functions import pheno_filt, img_filt, crop_img, norm_values, plot_pheno, plot_interactive
import cv2
import pandas as pd
import os
import pathlib
from loguru import logger
import re
from pandas.api.types import is_numeric_dtype

def img_match(data):
    
    print("\n################################# IMAGE MATCH ##################################\n")
    
    logger.info("Start image match process: This step will provide the plot of the phenotype(s) on the tiff image(s)\n")
     
    output=data["Paths"]["output_folder"]
            
    img_folder=os.path.join(data["Paths"]["data_input_folder"],"img_match") 
    csv_folder=os.path.join(data["Paths"]["data_input_folder"],"raw_data") 
    
    if len(data["Distance_match"]["pheno_list"])==0:
        pheno_list=data["Phenotypes"]["pheno_list"]
    else:
        pheno_list=data["Image_match"]["pheno_list"]

    output_f=os.path.join(output,"Img_match")
    pathlib.Path(output_f).mkdir(parents=True, exist_ok=True)
    
    dir=[f for f in os.listdir(img_folder) if not f.startswith('.')]
    img=list(set(list(map(lambda x: os.path.join(img_folder,x), dir))))
    csv=list(set(list(map(lambda x: os.path.join(csv_folder,x.split("_[")[0],re.search('.*?\]', x).group(0) +"_cell_seg_data.txt"), dir))))
    
    img.sort()
    csv.sort()
    
    for i,c in zip(img,csv):
        
        logger.info(f"Analyzing data {c}")
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
        
        #Filter by phenotype of interest
        pheno_df=pheno_filt(pheno_df,pheno_list)
        
        if len(pheno_df)==0:
            logger.warning(f"Data {c} has no match with requested Phenotype list. Skip to the next one!")
            continue
        
        #plot
        filename=os.path.splitext(os.path.basename(i))[0]+"_match_" + ''.join(map(str, pheno_list))
        plot_pheno(crop, pheno_df, os.path.join(output_f, filename))
        
        if data["Image_match"]["interactive"]:
            pathlib.Path(os.path.join(output_f, "Interactive_plots")).mkdir(parents=True, exist_ok=True) 
            plot_interactive(crop, pheno_df, os.path.join(os.path.join(output_f, "Interactive_plots"), filename), 
                             data["Image_match"]["layout_marker_edge_col"], 
                             data["Image_match"]["layout_marker_size"], 
                             data["Image_match"]["layout_xsize"], data["Image_match"]["layout_ysize"])
        
    logger.info("End image match!\n")
    
if __name__ == '__main__':
    img_match()