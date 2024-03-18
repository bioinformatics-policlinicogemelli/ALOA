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

#### Docker

3. Build docker file
```
cd <ALOA_folder_path>/ALOA
docker build -t aloa
```
4. Run docker in interactive mode
```
docker run -it -v <data_input_path>:./<data_input_folder_name> -v <output_data>:./output -v ./<path_to_config.json>:/config.json aloa
```

⚠️ *data_input_folder* inside json file and the mounted input folder must be the same. *i.e.* if data_input_folder="./data_input", the docker run command will be:
  ```
docker run -it -v <data_input_path>:./data_input -v <output_data>:./output -v ./<path_to_config.json>:/config.json aloa
```

## Usage
The first step to start using ALOA is to correctly set the configuration file *config.json*. This file is divided in 10 subsessions:
<br>
* **Paths**: here is possible to specify *data_input_folder*, *output_folder* and *sample_sheet* paths

* **Phenotypes**: here are specified the markers of interest into *pheno_list*

* **Descriptive**: here is possibile to specify parameters, for descriptive section, as *raw* and/or *normalized* count

* **Map_plot**: here are specified parameters as *multi_plot*, if you can... and *pheno_list* where specified the markers on performing the analysis

* **Distance**: here are specified the parameters for distance calculation as *save_csv* if you want to save the distance values, *save_img* uf you want to... and *pheno_list* if you want to calculate distance between specific markers

* **image_match**: here is possibile to specify a sublist of markers to print on composite images (if not specified the complete list of phenotypes will be considered). It is also possible to create interactive images, in addition to the static ones, by set the interactive option to true. Interactive graphs' layout can also be customize with the options layout_marker_edge_col and layout_marker_size to change respectively the edge color and the size of dots plotted on image and layout_xsize layout_ysize to customize the size of the image.

* **distance_match**: here is possibile to specify a sublist of markers whose distances are to be printed on composite images (if not specified the complete list of phenotypes will be considered).

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
All ALOA results are saved into a folder output whose path is specified in *config.json* file.

All log files, of this section, are saved into LOG folder with an intuitive name files.

For the merge, two types of files are genetared:
- files where **negative cellulas**, for markers of interest, haven't been deleted (saved into *Merged* folder)
- files where **negative cellulas**, for  markers of interest, have been deleted (saved into *Merged_clean* folder)

If there are two groups, into *Merged* and *Merged_clean*, will be created as many folders as groups present, containing merged files for each patient in the group. Each merged file is named as *patientnames.txt*

```
output_folder/
├── Log
|   ├── aloa_year_month_day_hour_minute_seconds.log
|   ├── clean_data_year_month_day_hour_minute_seconds.log
|   └── Merge_year_month_day_hour_minute_seconds.log
|
├── Merged
|   ├── Stroma
|   |   ├── Merge_cell_seg_data_Set4_1-6plex_S.txt
|   |   ├── ...
|   |    └──  Merge_cell_seg_data_Set12_20-6plex_S.txt
|   |
|   └──  Tumor
|       ├── Merge_cell_seg_data_Set4_1-6plex_T.txt
|       ├── ...
|       └──  Merge_cell_seg_data_Set12_20-6plex_T.txt
|
├── Merged_clean
|   ├── Stroma
|   |     ├── Merge_cell_seg_data_Set4_1-6plex_S.txt
|   |   ├── ...
|   |   └──  Merge_cell_seg_data_Set12_20-6plex_S.txt
|   |
|   └──  Tumor
|       ├── Merge_cell_seg_data_Set4_1-6plex_T.txt
|       ├── ...
|       └──  Merge_cell_seg_data_Set12_20-6plex_T.txt


```


### 2. Map Plots

### 3. Descriptive data + statistical analysis
From merged file, this section produces raw and/or normalized markers counts, for each patient and for each group, saved into csv files (*csv* folder) and visualize through **barplot** figures saved into *Barplot* folder.

⚠️ Two different formulas are used for normalized counts in barplots (*Norm_count_patientname.csv*) and in boxplots (*all_norm_count_patientsname*). For more details see [functions](./functions.md)


A *Descriptive* folder is created to save all the results of this section.

<p align="center"><img src="Image_readme/subfolder.png" width=320></p>


From raw/normalized counts, if there are two or more groups, a statistical comparison is made to understand if there are significance difference for markers count. The comparison is viewed through a **box plots** figure, with statistical annotation make though [TAP ](https://github.com/Discovery-Circle/tap)library, saved into *box_plot* folder

```
output_folder/
|    ├── Descriptive
|           ├── Stroma
|           |    ├── Bar_plot
|           |    |     ├── Bar_Plot_Normalized.jpeg
|           |    |     └── Bar_Plot_Raw.jpeg
|           |    └── csv
|           |          ├── all_norm_count_Set4_1-6plex_S.csv
|           |          ├── ...
|           |          ├── all_norm_count_Set12_20-6plex_S.csv
|           |          ├── Norm_count_Set4_1-6plex_S.csv
|           |          ├── ...
|           |          ├── Norm_count_Set12_20-6plex_S.csv
|           |          ├── Raw_count_Set4_1-6plex_S.csv
|           |          ├── ...
|           |          └── Raw_count_Set12_20-6plex_S.csv
|           |
|           ├── Tumor
|           |    ├── Bar_plot
|           |          ├── Bar_Plot_Normalized.jpeg
|           |          └── Bar_Plot_Raw.jpeg
|           |    ├── csv
|           |          ├── all_norm_count_Set4_1-6plex_T.csv
|           |          ├── ...
|           |          ├── all_norm_count_Set12_20-6plex_T.csv
|           |          ├── Norm_count_Set4_1-6plex_T.csv
|           |          ├── ...
|           |          ├── Norm_count_ Set12_20-6plex_T.csv
|           |          ├── Raw_count_Set4_1-6plex_T.csv
|           |          ├── ...
|           |          └──  Raw_count_Set12_20-6plex_T.csv                
|           |
|           |
|           └── Box Plots
|               ├── box_plot_comparison_Normalized.jpeg
|               └── box_plot_comparison_Rw.jpeg
                

```

### 4. Distances calculation

### 5. Satistical on distances

In this section a statistical analysis is performed to see if there are any significant differences, among the groups present, of the distributions of distances for the markers.
If *"pheno_from"* and *"pheno-to"* are specified in the configuration file, in , the analysis will be performed only for that type of distance. On the other hand, if **no** *"pheno_from"* and/or *"pheno-to"* is specified, the analysis is performed by permutation the different markers as "from" and "to".


A *Distance_Statistical* folder is created to save all the results of this section.

 <p align="center"><img src="Image_readme/distance_statistical.png" width=320> </p>

⚠️ Distance values are calculated as Z-score. For more details see For more details see [functions](./functions.md)



The results are: 
- **boxplot** figure,with statistical annotation make though [TAP ](https://github.com/Discovery-Circle/tap)library, saved in box_plot folder
- **csv** file with distance value for all patients for each group for a specific distance
- **distance_curve** figure (optional), saved into distance_curve folder, where are plotted the distributions of the distances of comparative groups for a specific couple of marker


```
output_folder/
|    ├── Distance_Statistical
|         ├── box_plot
|           |    ├── distances_box_plot_CD8+_to_CD68+.png
|           |    ├── ...
|           |    └── distances_box_plot_FoxP3+_to_CK+.png
|           ├── csv
|           |    ├── df_statistical_distance_CD8+_to_CD68+.csv
|           |    ├── ...
|           |    └── df_statistical_distance_FoxP3+_to_CK+.csv
|           |
|           ├── distance_curve
|           |    ├── plot_statistical_distance_CD8+_to_CD68+.png
|           |    ├── ...
|           |    └── plot_statistical_distance_FoxP3+_to_CK+.png
|           |
|           └── summary_statistical.csv

```

### 6. Clustering

### 7. Imaging


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
