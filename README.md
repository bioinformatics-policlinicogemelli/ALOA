# ALOA (Analysis spatiaL prOfiling imAging)

## Introduction

ALOA is a useful bioinformatics tool designed to transform raw data from PhenoCycler®-Fusion and inFormTM analysis ([Unlock the Power of Spatial Biology with phenoptrReports](https://www.akoyabio.com/wp-content/uploads/2022/01/Spatial-Biology-with-phentoprReports-TechNote.pdf)) into publication-ready results, thus advancing the accessibility and utility of spatial tissue analysis in cancer research.



## Features
- Spatial **Data** Analysis 📈

    This section provides an overall slide analysis and consists of several steps outlined below.

    The first step is the merge of cell seg data files of each roi. This files are merged into a single file for each patient. Cells that don't have any of the phenotypes of interest (OTHER) are removed. <br>

    The results will be saved into:<br>
    - <em><output_folder>/Merged_clean</em><br>



    From  *Merged_clean* files you can utilize 3 sections:
    - **Map Plot**: this section produces images to visualize spatial distribution of markers of interest
    
 
    - **Description**: this section procudes **bar plots**, through which visualize markers count for each patient, and **box plots** through which compare markers count from groups. <br> For the count there is the possibility to calculate: <br>
            - **RAW** count         <br>
            - **NORMALIZED** count <br>
    - **Distance**: this section calculated the distance from specific markers or from the combination of all the markers of the panel, starting from [phenoptr - Computing inter-cellular distances ](https://akoyabio.github.io/phenoptr/articles/computing_distances.html).<br> From distance data, it's possibile to make statistical analysis from groups, through box plot (default) and distribution curve plot (optional)

- Spatial **Imaging** Analysis <img src="./Image_readme/images.jpeg" width=20 height=12><br>

    This section provides a data visualization feature that produce static and interactive images from a selected number of ROIs. 

    Firstly ROIs' cell seg data are cleaned from cells that don't have any of the phenotypes of interest (OTHER).

    Subsequently, two different steps can be selected:
 
    - **Image Match**: This section provides a spatial plot of phenotypes' position on ROI's composite image.
    The results will be saved into <em><output_folder>/Img_match</em><br>

    - **Distance Match**: This section provides a spatial plot of phenotypes' distances on ROI's composite image for all possible combination of couples of the selected phenotypes. The results will be saved into <em><output_folder>/Distance_match</em><br>



## Installation


1. Open a terminal
2. Digit the following command to clone the repository folder: 
```
git clone https://github.com/bioinformatics-policlinicogemelli/ALOA
```
3.  Install all of the packages required
```
cd <ALOA_folder_path>/ALOA

pip install -r requirements.txt

Rscript installation_rpackages.R req.txt
```

From Docker

## Usage
The first step to start using ALOA is to correctly set the configuration file *config.json*. This file is divided in 10 subsessions:
<br>
* **Paths**: here is possible to specify *data_input_folder*, *output_folder* and *sample_sheet* paths

* **Packages**: 

* **Phenotypes**: here are specified the markers of interest into *pheno_list*

* **Descriptive**: here is possibile to specify parameters, for descriptive section, as *raw* and/or *normalized* count

* **Map_plot**: here are specified parameters as *multi_plot*, if you can... and *pheno_list* where specified the markers on performing the analysis

* **Distance**: here are specified the parameters for distance calculation as *save_csv* if you want to save the distance values, *save_img* uf you want to... and *pheno_list* if you want to calculate distance between specific markers

* **image_match**:

* **distance_match**:

* **statistical_distance**: here is possibile to specify the markers for which you can perform the distance statsical analysis in *pheno_from* and *pheno_to*

* **cluster**: here is possibile to specify the parameters for clustering analysis as *pheno_list* if you want to specify specific markers, *k* if you want to specify the number of clusters, *cluster_method* to choose the clustering method (spectral or kmeans)

#INSERIRE WORFLOW
### 1. Merge cell seg data
In this section the cell_seg_data.txt files, for each patient, are merged into a single file.
The single cell_seg_data are in patient specific folder into data_test/raw_data.
```
data_test/
├── img_match
│   ├── Set4_1-6plex_[11472,51360]_cell_seg_data.txt
│   ├── Set4_1-6plex_[11472,51360]_composite_image.jpeg
|   ├── ...
|   ├── Set12_20-6plex_[17241,54367]_cell_seg_data.txt
│   └── Set12_20-6plex_[17241,54367]_composite_image.jpeg
|
├── raw_data
│   ├── Set4_1-6plex_S
|       ├── Set4_1-6plex_[11472,51360]_cell_seg_data.txt
|       ├── ...
|       └──  Set4_1-6plex_[16142,55840]_cell_seg_data.txt
|   ├── ...
│   └── Set12_20-6plex_S
|       ├── Set12_20-6plex_[14146,53503]_cell_seg_data.txt
|       ├── ...
|       └──  Set12_20-6plex_[17241,54367]_cell_seg_data.txt

|
└──  sample_sheet.xlsx

```
The merged files are saved in a specific folder (*Merged*) into the output folder specified in the configuration file. If there are two groups, into *Merged* will be created as many folders as groups present, containing merged files for each patient in the group. Each merged file is named as *patientnames.txt*

```
output_folder/
├── Merged
|   ├── Stroma
|       ├── Set4_1-6plex_S.txt
|       ├── ...
|       └──  Set12_20-6plex_S.txt
|
|   └──  Tumor
|       ├── Set4_1-6plex_T.txt
|       ├── ...
|       └──  Set12_20-6plex_T.txt
|
```


### 2. Map Plots

### 3. Descriptive data + statistical analysis
From merged file, this section produces raw and/or normalized markers counts, for each patient and for each group, saved into csv files and visualize through a **barplot**

A *Descriptive* folder is created to save all the results of this section. For major details see [functions](./functions.md)

<p align="center"><img src="Image_readme/subfolder.png" width=320></p>




From raw/normalized counts, if there are two or more groups, a statistical comparison is made to understand if there are significance difference for markers count. The comparison is viewed through a **box plots** with statistical annotation make though [TAP ](https://github.com/Discovery-Circle/tap)library, saved into *box_plot* folder

```
output_folder/
|    ├── Descriptive
|           ├── Group Stroma
|           |    ├── Barplot
|           |    |     ├── Bar_Plot_Raw.jpeg
|           |    |     └──  Bar_Plot_Normalized.jpeg
|           |    └── csv
|           |          ├── Norm_count_Set4_1-6plex_S.csv
|           |          ├── ...
|           |          ├── Norm_count_ Set12_20-6plex_S.csv
|           |          ├── Rwe_count_Set4_1-6plex_S.csv
|           |          ├── ...
|           |          └──  Raw_count_Set12_20-6plex_S.csv
|           |
|           ├── Group Tumor
|           |    ├── Barplot
|           |          ├── Bar_Plot_Raw.jpeg
|           |          └──  Bar_Plot_Normalized.jpeg
|           |    ├── csv
|           |          ├── Norm_count_Set4_1-6plex_T.csv
|           |          ├── ...
|           |          ├── Norm_count_ Set12_20-6plex_T.csv
|           |          ├── Raw_count_Set4_1-6plex_T.csv
|           |          ├── ...
|           |          └──  Raw_count_Set12_20-6plex_T.csv                
|           |
|           |
|           └──  Box_plot

```

### 4. Distances calculation + statistical analysis

### 5. Imaging


## Launch ALOA main

These are the options that can be set for this block:

| Options | Input | Type | Required
|----------------|----------------| :---:| :---:|
|-m <br> --merge| <p align="justify">merge datasets from ROIs of the same patient| boolean | No
|-M <br> --maps| <p align="justify">create maps plot| boolean | No
|-d <br> --distance| <p align="justify">distance evaluation between phenotypes| boolean | No
|-s <br> --stats| <p align="justify">create distance stats| boolean | No
|-o <br> --overview| <p align="justify">create data overview | boolean | No
|-C <br> --clustering| <p align="justify">cluster data| boolean | No
|-I <br> --clustering| <p align="justify">create phenotypes image match| boolean | No
|-D <br> --clustering| <p align="justify">create phenotypes distance match| boolean | No
|-w <br> --clustering| <p align="justify">create data overview| boolean | No
