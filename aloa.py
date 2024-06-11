import os
import json
import argparse
from loguru import logger
import logging
import pathlib
import sys
import rpy2.robjects as ro
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
rpy2_logger.setLevel(logging.ERROR)
from datetime import datetime
import glob
import pandas as pd

from function_descriptive_analysis import main as descriptive
from functions_statistical_distance import main as stat_dist
from clustering import main as clust
from img_match import img_match
from distance_match import distance_match
from cross_PCF import main as pcf

@logger.catch()

#*****************************************************************

def log_settings(logout):
    logfile="aloa_"+datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+".log"
    logger.level("INFO", color="<green>")
    logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True, catch=True, backtrace=True, diagnose=True)
    logger.add(os.path.join(logout,logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w", backtrace=True, diagnose=True)
    
    return(os.path.join(logout,logfile))
    
#*****************************************************************

def all_true(args):
    for arg in vars(args):
        if arg=="all":
            continue
        print(arg, getattr(args, arg))
        setattr(args, arg,True)
        print(arg, getattr(args, arg))

#*****************************************************************
        
def check_output(output,fun):
    if fun=="merge" or fun=="dist":
        path_to_find=os.path.join(output,"*","*txt")
    elif fun=="maps":
        path_to_find=os.path.join(output,"*","*pdf")
    else:
        return("")
    check=glob.glob(path_to_find)
    return(check)

#*****************************************************************

def merge_log(log_list):

    os.system("cat " + log_list[1] + " >> " + log_list[0])
    os.system("rm " + log_list[1])
    logger.remove()
    logger.add(log_list[0],format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="a", backtrace=True, diagnose=True)        
    
#*****************************************************************

def check_log(log, n_warn=0, n_err=0, n_crit=0):
    
    log_txt=open(log, 'r').read()
    n_warn+=int(log_txt.count("WARNING"))
    n_err+= int(log_txt.count("ERROR"))
    n_crit+=int(log_txt.count("CRITICAL"))
    return(n_warn,n_err,n_crit)

#*****************************************************************
#*****************************************************************

def aloa(args, data, logfile):
    
    logger.info(f"aloa.py args: [merge:{args.merge}, mapping:{args.maps}, distance:{args.distance}, img_match:{args.imgMatch}, dst_match:{args.dstMatch}, overview:{args.overview}, stats:{args.stats}, cluster:{args.clustering}, all:{args.all}]")
    
    output=data["Paths"]["output_folder"]
    
    if not(len(os.listdir(output))==1 and os.listdir(output)[0]=="Log"):
        logger.critical(f"It seems that {output} folder already exists! Delete the folder or change the output name in the config file!")
        return()
        
    n_err=0
    n_warn=0
    n_crit=0
             
    #########################################
    #          DATA REORGANIZATION          #
    #########################################
    
    if args.merge:
        
        logger.info("|-> Merge step starting now")
        ro.r['source']('merge.R')
        merge = ro.globalenv['merge']
        log_name_merge=merge()[0]
        merge_log([logfile,log_name_merge])
        
        logger.info("|-> Clean step starting now")
        ro.r['source']('clean_data.R')
        clean = ro.globalenv['clean']
        log_name_clean=clean()[0]
        merge_log([logfile,log_name_clean])
                           
        #########################################
        #             DESCRIPTIVE               #
        #########################################
        
        if args.overview:
            logger.info("|-> Descriptive step starting now\n")
            descriptive()
        
        #########################################
        #             MAPPING PHENO             #
        ######################################### 
        
        if args.maps:
            logger.info("|-> Maps plot step starting now")
            ro.r['source']('maps_plot.R')
            maps = ro.globalenv['maps']
            log_name=maps()[0]
            merge_log([logfile,log_name])
            
        #########################################
        #               DISTANCE                #
        #########################################
        
        if args.distance:
            logger.info("|-> Distance evaluation step starting now")
            ro.r['source']('distance_eval.R')
            distance = ro.globalenv['distance']
            log_name=distance()[0]
            merge_log([logfile,log_name])
            
            if args.stats:
                logger.info("|-> Distance statistical evaluation step starting now")
                stat_dist()
                
        #########################################
        #              CLUSTERING               #
        ######################################### 
        
        if args.clustering:
            logger.info("|-> Clustering step starting now")
            clust()
           
    #########################################
    #               IMAGING                 #
    #########################################        
    if args.imgMatch:          
        logger.info("|-> Image matching step starting now")
        img_match()
                       
    if args.dstMatch: 
        logger.info("|-> Distance matching step starting now")
        distance_match()
    
    #########################################
    #                  PCF                  #
    #########################################

    if args.pcf:
        logger.info("|-> PCF step starting now")
        pcf()
    
    n_warn, n_err, n_crit=check_log(logfile)
    logger.info(f"ALOA script completed with {n_warn} warnings, {n_err} errors and {n_crit} critical!")
    
#### 
class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)
#### 

def main(): 
    
    with open("config.json") as f:
        data=json.load(f)
        
    output=data["Paths"]["output_folder"]
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    log_path=os.path.join(output,"Log")
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True) 

    logger.remove()
    
    logfile=log_settings(log_path)
    
    logger.info("Welcome to ALOA")
         
    parser = MyArgumentParser(add_help=True, usage=None, description='Argument of ALOA script')
    
    # MERGE + CLEAN BLOCK
    parser.add_argument('-m', '--merge', required=False, action='store_true', help='merge datasets from ROIs of the same patient')
    
    # MAPS PLOT BLOCK
    parser.add_argument('-M', '--maps', required=False, action='store_true', help='create maps plot')
    
    # DISTANCE BLOCK
    parser.add_argument('-d', '--distance', required=False, action='store_true', help='distance evaluation between phenotypes')
    parser.add_argument('-s', '--stats', required=False, action='store_true', help='create distance stats')
    
    # STATISTICAL BLOCK
    parser.add_argument('-o', '--overview', required=False, action='store_true', help='create data overview')
    
    # CLUSTERING BLOCK
    parser.add_argument('-c', '--clustering', required=False, action='store_true', help='cluster data')

    # PCF BLOCK
    parser.add_argument('-p', '--pcf', required=False, action='store_true', help='cluster data')
        
    # IMAGING BLOCK
    parser.add_argument('-I', '--imgMatch', required=False, action='store_true', help='create phenotypes image match')
    parser.add_argument('-D', '--dstMatch', required=False, action='store_true', help='create phenotypes distance match')
    
    # ALL
    parser.add_argument('-a', '--all', required=False, action='store_true', help='do all analysis')
    
    try:
        args = parser.parse_args()  
    except ValueError as e:
        logger.critical(f"Argument error: {e}!")
        exit(1)
    args.merge=True
    args.maps=True
    args.distance=True     
    if args.all:
        all_true(args)
        
    if (not args.imgMatch and not args.dstMatch and not args.pcf) and not args.merge:
        logger.critical("This pipeline requires the merge (-m) option to proceed with the total analysis. If not setted, at least one of the single image options (-I, -D or -p) must be select. Check your input options")
        exit(1)
    
    if args.stats and not args.distance:
        logger.critical(f"It seems that stats (-s) option is provided without the distance (-d) one. That is not possible beacause statistical evaluation requires distance evaluation. Try launch the command again setting -d and -s concurrently.")
        exit(1)
            
    aloa(args, data, logfile)

if __name__ == '__main__':
    main()
    