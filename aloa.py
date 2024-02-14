import os
import json
import argparse
import logging
import rpy2.robjects as ro
from function_descriptive_analysis import main as descriptive
from functions_statistical_distance import main as stat_dist
import pathlib

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger('aloa')

def aloa(args):
    
    f=open("config.json")
    data=json.load(f)
    output=data["Paths"]["output_folder"]
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    log_path=os.path.join(output,"Log")
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True) 
    fh = logging.FileHandler(os.path.join(log_path,'aloa.log'))
    # ch = logging.StreamHandler()
    # ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    # ch.setFormatter(formatter)
    # logging.getLogger().setLevel(logging.INFO)

    logger.info(f"aloa.py args [merge:{args.merge}, clear:{args.clean},\
                mapping:{args.maps},distance:{args.distance},\
                img_match:{args.imgMatch}, dst_match:{args.dstMatch},\
                overview:{args.overview}, stats:{args.stats},\
                cluster:{args.clustering},\
                overwrite:{args.overwrite}]")
    
    if args.overwrite and os.path.isdir(output):
        logger.info(f"{output} folder already exists and will be overwritten")
        
    #########################################
    #          DATA REORGANIZATION          #
    #########################################
                    
    if args.merge:
        logger.info("Merge step starting now")
        ro.r.source("./merge.R", encoding="utf-8")
        
        if args.clean:
            if args.merge or os.listdir(os.path.join(data["Paths"]["output_folder"],"Merged")>0):
                logger.info("Clean step starting now")
                ro.r.source("./clean_data.R", encoding="utf-8")
            else:
                logger.error("Merge folder empty")
        
        #########################################
        #             DESCRIPTIVE               #
        ######################################### 
        
        if args.overview:
            logger.info("Descriptive step starting now")
            descriptive()
            
        #########################################
        #             MAPPING PHENO             #
        #########################################  
        
        if args.maps:
            logger.info("Maps plot step starting now")
            ro.r.source("./maps_plot.R", encoding="utf-8")

        #########################################
        #               DISTANCE                #
        #########################################  

        if args.distance:
            logger.info("Distance evaluation step starting now")
            ro.r.source("./distance_eval.R", encoding="utf-8")
            
            if args.stats:
                logger.info("Distance statistical evaluation step starting now")
                stat_dist()
                
        #########################################
        #             CLUSTERING                #
        ######################################### 
        
        if args.clustering:
            logger.info("Clustering step starting now")

    #########################################
    #               IMAGING                 #
    #########################################        
    if args.imgMatch:          
        logger.info("Image matching step starting now")
         
    if args.dstMatch: 
        logger.info("Distance matching step starting now") 
    
#### 
class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

def main(): 
     
    parser = MyArgumentParser(add_help=True, exit_on_error=True, usage=None, description='Argument of ALOA script')
    
    # MERGE + CLEAN BLOCK
    parser.add_argument('-m', '--merge', required=False, action='store_true', help='merge datasets from ROIs of the same patient')
    parser.add_argument('-c', '--clean', required=False, action='store_true', help='remove the all Phenotype-negative cells')
    
    # MAPS PLOT BLOCK
    parser.add_argument('-M', '--maps', required=False, action='store_true', help='create maps plot')
    
    # DISTANCE BLOCK
    parser.add_argument('-d', '--distance', required=False, action='store_true', help='distance evaluation between phenotypes')
    parser.add_argument('-s', '--stats', required=False, action='store_true', help='create distance stats')
    
    # STATISTICAL BLOCK
    parser.add_argument('-o', '--overview', required=False, action='store_true', help='create data overview')
    
    # CLUSTERING BLOCK
    parser.add_argument('-C', '--clustering', required=False, action='store_true', help='cluster data')
        
    # IMAGING BLOCK
    parser.add_argument('-I', '--imgMatch', required=False, action='store_true', help='create phenotypes image match')
    parser.add_argument('-D', '--dstMatch', required=False, action='store_true', help='create phenotypes distance match')
    
    # overwrite folder
    parser.add_argument('-w', '--overwrite', required=False, action='store_true', help='create data overview')
    
    try:
        args = parser.parse_args()  
    except ValueError :
        logger.critical("Something went wrong! Check your inputs")
        exit(1)
        
    if (not args.imgMatch or not args.dstMatch) and not args.merge:
        logger.critical("Something went wrong! This pipeline requires the merge (-m) option to proceed with the total analysis or one of the single image options (-I or -D)")
        exit(1)
    
    aloa(args)
    
if __name__ == '__main__':
    main()