import os
import json
import argparse
import logging
import pathlib
from datetime import datetime
import rpy2.robjects as ro
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
rpy2_logger.setLevel(logging.ERROR)

from function_descriptive_analysis import main as descriptive
from functions_statistical_distance import main as stat_dist
from clustering import main as clust
from img_match import img_match
from distance_match import distance_match

def aloa(args, data):
    
    logging.info(f"\naloa.py args: \n[merge:{args.merge}, clear:{args.clean},\
                 \nmapping:{args.maps},distance:{args.distance},\
                 \nimg_match:{args.imgMatch}, dst_match:{args.dstMatch},\
                 \noverview:{args.overview}, stats:{args.stats},\
                 \ncluster:{args.clustering},\
                 \noverwrite:{args.overwrite}]")
    
    output=data["Paths"]["output_folder"]
    
    if args.overwrite and os.path.isdir(output):
        logging.info(f"{output} folder already exists and will be overwritten")
    
    #########################################
    #          DATA REORGANIZATION          #
    #########################################
    
    if args.merge:
        logging.info("Merge step starting now")
        ro.r.source("./merge.R", encoding="utf-8")
        
        if args.clean:
            if args.merge or os.listdir(os.path.join(data["Paths"]["output_folder"],"Merged")>0):
                logging.info("Clean step starting now")
                ro.r.source("./clean_data.R", encoding="utf-8")
            else:
                logging.error("Merge folder empty")

        if not args.clean:
            logging.warning("Clean step was not selected. This step is strongly recommended")
        
        #########################################
        #             DESCRIPTIVE               #
        #########################################
        
        if args.overview:
            logging.info("Descriptive step starting now")
            descriptive()
        
        #########################################
        #             MAPPING PHENO             #
        #########################################  
        
        if args.maps:
            logging.info("Maps plot step starting now")
            ro.r.source("./maps_plot.R")
            
        #########################################
        #               DISTANCE                #
        #########################################  

        if args.stats and not args.distance:
            logging.info("It seems that stats (-s) option is provided without the distance (-d) one. \
                        That is not possible beacause statistical evaluation require distance evaluation. \
                            Try launch the command again setting -d and -s concurrently.")
            exit()
        
        if args.distance:
            logging.info("Distance evaluation step starting now")
            ro.r.source("./distance_eval.R", encoding="utf-8")
            
            if args.stats:
                logging.info("Distance statistical evaluation step starting now")
                stat_dist()
                
        #########################################
        #             CLUSTERING                #
        ######################################### 
        
        if args.clustering:
            logging.info("Clustering step starting now")
            clust()
            
    #########################################
    #               IMAGING                 #
    #########################################        
    if args.imgMatch:          
        logging.info("Image matching step starting now")
        img_match()
         
    if args.dstMatch: 
        logging.info("Distance matching step starting now")
        distance_match()
    
#### 
class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

def main(): 
    
    f=open("config.json")
    data=json.load(f)
    output=data["Paths"]["output_folder"]
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    log_path=os.path.join(output,"Log")
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True) 

    format_time=datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    logging.basicConfig(format="[%(levelname)s] %(asctime)s - %(message)s", datefmt='%Y-%m-%d %H:%M:%S',filename=os.path.join(log_path,"aloa_"+format_time+".log"),filemode="a")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    logger = logging.getLogger()  # get the root logger
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("[%(levelname)s] %(asctime)s - %(message)s", "%Y-%m-%d %H:%M:%S"))
    logger.addHandler(console_handler)
         
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
    
    if (args.clean and not args.merge):
        logger.critical("Something went wrong! This pipeline requires the merge (-m) option before the clean (-c) one!")
        exit(1)
       
    if (not args.imgMatch or not args.dstMatch) and not args.merge:
        logger.critical("Something went wrong! This pipeline requires the merge (-m) option to proceed with the total analysis or one of the single image options (-I or -D)")
        exit(1)
    
    aloa(args, data)
    
if __name__ == '__main__':
    main()