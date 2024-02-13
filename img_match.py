#libraries import
from image_proc_functions import pheno_filt, img_filt, crop_img, norm_values, plot_pheno, plot_interactive
import cv2
import pandas as pd
import json
import os
import pathlib

def img_match():
    with open('config.json', 'r') as myfile:
        data=myfile.read()
    # parse json file
    obj = json.loads(data)

    output=obj["Paths"]["output_folder"]
    pathlib.Path(output).mkdir(parents=True, exist_ok=True) 
    output_f=os.path.join(output,"Img_match")
    pathlib.Path(output_f).mkdir(parents=True, exist_ok=True)

    input_folder=os.path.join(obj["Paths"]["data_input_folder"],"img_match") 
    #pheno_list=obj["image_match"]["pheno_list"]
    if len(obj["distance_match"]["pheno_list"])==0:
        pheno_list=obj["Phenotypes"]["pheno_list"]
    else:
        pheno_list=obj["image_match"]["pheno_list"]

    dir=[f for f in os.listdir(input_folder) if not f.startswith('.')]
    dir=list(set(list(map(lambda x: os.path.join(input_folder,x.split(']_')[0]+"]"), dir))))
    
    for d in dir:
        
        print("Analyzing data "+ d)
        
        #dataframe section
        df=pd.read_csv(d+"_cell_seg_data.txt",sep="\t")
        
        # image section
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
        
        if obj["image_match"]["interactive"]:
            pathlib.Path(os.path.join(output_f, "Interactive_plots")).mkdir(parents=True, exist_ok=True) 
            plot_interactive(crop, pheno_df, os.path.join(os.path.join(output_f, "Interactive_plots"), filename), 
                             obj["image_match"]["layout_marker_edge_col"], 
                             obj["image_match"]["layout_marker_size"], 
                             obj["image_match"]["layout_xsize"], obj["image_match"]["layout_ysize"])
        
    print("All done!")
    
if __name__ == '__main__':
    img_match()