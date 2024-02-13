#libraries import
import os
import numpy as np
import cv2
from matplotlib import pyplot as plt
import pandas as pd
pd.DataFrame.iteritems = pd.DataFrame.items

import seaborn as sns

from skimage.metrics import structural_similarity

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

import itertools

from chart_studio import plotly as py
import plotly.graph_objs as go
import plotly.express as px
from PIL import Image

### Image Preprocessing Section ###

# function to filter dataframe (remove phenotype not of interest)
def pheno_filt(pheno_df, pheno_list=[]):    
    pheno_col=pheno_df[pheno_df.columns[pd.Series(pheno_df.columns).str.startswith('Phenotype')]]
    
    pheno_col=pheno_col.columns
    pheno_df.fillna('other', inplace=True)
    pheno_df.fillna("OTHER", inplace=True)
    pheno_df=pheno_df.apply(lambda x: x.replace("other",""))
    pheno_df=pheno_df.apply(lambda x: x.replace("OTHER",""))
    pheno_df["Pheno"]=pheno_df[pheno_col].sum(axis=1)
    pheno_df=pheno_df[pheno_df['Pheno']!=""]
    
    if not pheno_list:
        return pheno_df
    
    pheno_list_mod=[]
    for ph in set(pheno_list):
        if len(set(pheno_df["Pheno"].values).intersection(set([ph])))>0:
            pheno_list_mod.append(list(set(pheno_df["Pheno"].values).intersection(set([ph])))[0])
        elif len(ph.split(','))>1:
            combinations=list(itertools.permutations(ph.split(',')))
            for cc in combinations:
                cc_join=''.join(cc)
                if len(set(pheno_df["Pheno"].values).intersection(set([cc_join])))>0:
                    pheno_list_mod.append(cc_join)
    
    idx=[]
    for cc in pheno_list_mod:
        idx.append(pheno_df.index[pheno_df['Pheno'] == cc].tolist())
    idx=sum(idx, [])
    pheno_df=pheno_df.loc[idx,:]
    
    return pheno_df.reset_index()

# function to filter image (remove isolated pixels)
def img_filt(img):
    blur = cv2.GaussianBlur(img, (15, 15), 2)
    hsv = cv2.cvtColor(blur, cv2.COLOR_BGR2HSV)
    lower_green = np.array([37, 0, 0])
    upper_green = np.array([179, 255, 255]) #179
    mask = cv2.inRange(hsv, lower_green, upper_green)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (15, 15))
    opened_mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
    masked_img = cv2.bitwise_and(img, img, mask=opened_mask)
    
    #diff=img-masked_img
    img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    masked_img_gray = cv2.cvtColor(masked_img, cv2.COLOR_BGR2GRAY)
    
    #evaluate image similarity between masked and figure. If similarity is low,
    #filter has masked too much so the original image is returned
    score, _ = structural_similarity(img_gray, masked_img_gray, full=True)
    
    if score<.4:
        print("Low similiraty: {:.2f}%".format(score))
        return img
    
    return masked_img

# function to crop images
def crop_img(image,th):
    y_nonzero, x_nonzero, _ = np.nonzero(image>th)
    return image[np.min(y_nonzero):np.max(y_nonzero), np.min(x_nonzero):np.max(x_nonzero)]

#function to normalize X and Y cell positions
def norm_values(pheno_df, crop):
    pheno_df["Cell X Position Norm"]=(
        (pheno_df["Cell X Position"]-min(pheno_df["Cell X Position"]))/
        (max(pheno_df["Cell X Position"])-min(pheno_df["Cell X Position"]))
        )*crop.shape[1]

    pheno_df["Cell Y Position Norm"]=(
        (pheno_df["Cell Y Position"]-min(pheno_df["Cell Y Position"]))/
        (max(pheno_df["Cell Y Position"])-min(pheno_df["Cell Y Position"]))
        )*crop.shape[0]

    return pheno_df

#function to plot pheno images
def plot_pheno(img, pheno_df, filename , p='gist_ncar', edge_col='black', dim_dot=15):
    implot=plt.imshow(img)
    sns.scatterplot(pheno_df, 
                    x=pheno_df["Cell X Position Norm"], 
                    y=pheno_df["Cell Y Position Norm"],
                    s=dim_dot,
                    hue="Pheno",
                    palette=p,
                    edgecolor=edge_col)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, markerscale=2)
    plt.xlabel("X Cell Position")
    plt.ylabel("Y Cell Position")
    
    plt.savefig(filename+".tif",dpi=300, format="tiff", bbox_inches='tight')
    
    plt.close()

### Distance Section ###

#R environment settings
r = ro.r
r['source']('distance_match_func.R')
distance_eval_fun_r = ro.globalenv['distance_eval']

#function to evaluate distance    
def dist_eval(pheno_df):
    with localconverter(ro.default_converter + pandas2ri.converter):
            #convert python df to R df
        pheno_df_r = ro.conversion.py2rpy(pheno_df)
            #launch r function
    df_result_r = distance_eval_fun_r(pheno_df_r)
            #convert R df to python df
    pheno_df=ro.pandas2ri.rpy2py_dataframe(df_result_r)
    return pheno_df

#function to plot distance images
def plot_dist(img, pheno_df, ph1_to_ph2, filename, dim_dot=8, m=False):
    
    ph1_to_ph2=ph1_to_ph2.reset_index()
    
    ph1=ph1_to_ph2["Phenotype"][0]
    ph2=[col for col in ph1_to_ph2 if col.startswith('Phenotype.')][0].split('.')[-1]
    
    ph1_df=pheno_df[pheno_df["Phenotype"]==ph1].reset_index()
    ph2_df=pheno_df[pheno_df["Phenotype"]==ph2].reset_index()
        
    implot=plt.imshow(img)
    
    plt.plot(ph1_df['Cell X Position Norm'], ph1_df['Cell Y Position Norm'],'k.', markerfacecolor='c', markersize=dim_dot, label=ph1)
    plt.plot(ph2_df['Cell X Position Norm'], ph2_df['Cell Y Position Norm'],'k.', markerfacecolor='r', markersize=dim_dot, label=ph2)
    
    for n in range(0, len(ph1_to_ph2), 2):
        plt.plot([ph1_to_ph2.loc[n, 'Cell X Position Norm'], ph1_to_ph2.loc[n, 'Cell X Position Norm.'+ph2]], 
                 [ph1_to_ph2.loc[n, 'Cell Y Position Norm'], ph1_to_ph2.loc[n, 'Cell Y Position Norm.'+ph2]], 'w--', alpha=0.8)
                
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, markerscale=2)
    plt.xlabel("X Cell Position ($\mu$m)")
    plt.ylabel("Y Cell Position ($\mu$m)")
    
    if not m:
        plt.title("Nearest "+ph1+" to each "+ph2)
        plt.savefig(filename+"_Nearest_"+ph1+"_to_each_"+ph2+".tif",dpi=300, format="tiff", bbox_inches='tight')
    else:
        plt.title("Mutual nearest neighbors -"+ph1+" and "+ph2)
        plt.savefig(filename+"_Mutual_nearest_neighbors_"+ph1+"_and_"+ph2+".tif",dpi=300, format="tiff", bbox_inches='tight')
    
    plt.close()

def plot_interactive(img, pheno_df, filename , edge_col='black', dim_dot=10, xsize_img=1200, ysize_img=900):
    image=Image.fromarray(img)

    data=[]
    for group in pheno_df["Pheno"].unique():
        df_group=pheno_df[pheno_df["Pheno"]==group]
        trace=go.Scatter(x=df_group["Cell X Position Norm"], y=-df_group["Cell Y Position Norm"],
                         mode='markers',
                         marker=dict(size=dim_dot,line=dict(width=1,color=edge_col)),
                         name=group)
        data.append(trace)
        
    layout=go.Layout(title="Phenotype(s)", width=xsize_img, height=ysize_img)
    
    fig=go.Figure(data=data,layout=layout)
    
    fig.add_layout_image(
                  source= image,
                  xref= "x",
                  yref= "y",
                  x= 0,
                  y= 0,
                  xanchor='left',
                  yanchor='top',
                  sizex= image.width,
                  sizey= image.height,
                  sizing= "stretch",
                  opacity= 1,
                  layer= "below")
    
    fig.update_layout(xaxis=dict(showgrid=False), 
                      yaxis=dict(showgrid=False), 
                      plot_bgcolor='white', 
                      xaxis_range=[0, image.width], yaxis_range=[-image.height, 0],
                      xaxis_title="", yaxis_title="",
                      xaxis_showticklabels=False, yaxis_showticklabels=False,
                      showlegend=True
                      )
    #fig.show()
    fig.write_html(filename+".html")