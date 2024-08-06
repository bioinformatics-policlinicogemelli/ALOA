# ALOA (Analysis spatiaL prOfiling imAging)

## Introduction

ALOA is a useful bioinformatics tool designed to transform raw data from PhenoCycler®-Fusion and inFormTM analysis ([Unlock the Power of Spatial Biology with phenoptrReports](https://www.akoyabio.com/wp-content/uploads/2022/01/Spatial-Biology-with-phentoprReports-TechNote.pdf)) into publication-ready results, thus advancing the accessibility and utility of spatial tissue analysis in cancer research.



## Features
- Spatial **Data** Analysis 📈

    This section provides an overall slide analysis and consists of several steps outlined below.

    The first step is the merge of cell seg data files of each roi. These files are merged into a single file for each patient. Cells that don't have any of the phenotypes of interest (OTHER) are removed. <br>

    The results are saved into:<br>
    - <em><output_folder>/Merged_clean</em><br>



    Starting from  *Merged_clean* files it is possible to proceed for 5 different sections:
    - **Map Plot**: this section produces images to visualize spatial distribution of markers of interest

    - **Description**: this section produces **bar plots**, through which visualize markers row and normalized count for each patient, and **box plots** through which compare markers count from groups. <br> For the count there is the possibility to calculate: <br>

    - **Distance**: this section integrates [find_nearest_distance](https://akoyabio.github.io/phenoptr/articles/computing_distances.html) phenoptr script and calculate the distance between couples of markers.<br> Distances can be compared between groups and visualize through box plot (default) and distribution curve plot (optional)

    - **Clustering**: this section evaluate different clustering algorithm on data. Markers distribution across clusters can be seen and analyze through boxplot and reports

    - **CrossPCF**: this section integrates [crossPCF](https://github.com/JABull1066/ExtendedCorrelationFunctions) script from Bull et al, 2024 and calculate the cross-PCF between couples of markers. <br> A statistical analysis can be conducted if two or more groups are specified.

- Spatial **Imaging** Analysis <img src="./Image_readme/images.jpeg" width=20 height=12><br>

    This section provides a data visualization feature that produce static and interactive images from a selected number of ROIs. 

    Firstly ROIs' cell seg data are cleaned from cells that don't have any of the phenotypes of interest (OTHER).

    Subsequently, 3 different options can be selected:
 
    - **Image Match**: This section provides a spatial plot of phenotypes' position on ROI's composite image.
    
    - **Distance Match**: This section provides a spatial plot of phenotypes' distances on ROI's composite image for all possible combination of couples of the selected phenotypes.
    
    - **Cross-PCF Image Match**: This section provides a termic plot of Cross-PCF on ROI's composite image for all possible combination of couples of the selected phenotypes.


## Installation and usage

### Docker (recommended)

#### Requirements:
[Docker ](https://www.docker.com/)
#### Procedure:

1. Open a terminal
2. Clone the repository folder:
```
git clone https://github.com/bioinformatics-policlinicogemelli/ALOA
```
3. Build docker file
```
cd <ALOA_folder_path>/ALOA
docker build -t aloa .
```
4. Test installation
```
docker run -it -v <input_folder>:/input -v <output_folder>:/output aloa
:/# python aloa.py -a 
```

### Local

#### Requirements:
[Python ](https://www.python.org/)>= 3.10 <br>
[R ](https://www.r-project.org/)>= 4.3 <br>
⚠️ for Windows users [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installation is also required

#### Procedure:

1. Open a terminal
2. Clone the repository folder:
```
git clone https://github.com/bioinformatics-policlinicogemelli/ALOA
```
3.  Install all of the packages required
```
cd <ALOA_folder_path>/ALOA

pip install -r requirements.txt

Rscript installation_rpackages.R req.txt
```

4. Test installation
```
python aloa.py -a
```

## Options

These are the options that can be set by user:

| Options | Input | Type | Required
|----------------|----------------| :---:| :---:|
|-m <br> --merge| <p align="justify">merge datasets from ROIs of the same patient| boolean | Yes*
|-M <br> --maps| <p align="justify">create maps plot| boolean | No
|-d <br> --distance| <p align="justify">distance evaluation between phenotypes| boolean | No
|-s <br> --stats| <p align="justify">create distance stats| boolean | No
|-o <br> --overview| <p align="justify">create data overview | boolean | No
|-c <br> --clustering| <p align="justify">cluster data| boolean | No
|-p <br> --pcf| <p align="justify">pcf analysis| boolean | No
|-I <br> --imgMatch| <p align="justify">create phenotypes image match| boolean | No
|-D <br> --dstMatch| <p align="justify">create phenotypes distance match| boolean | No
|-a <br> --all| <p align="justify">do all the analysis| boolean | No

*only for whole slide analysis (**NOT** required for pcf, imgMatch and dstMatch)


## Usage
The first step to start using ALOA is to correctly set the configuration file *config.json*. This file is divided in 11 subsessions:

* **Paths**: here is possible to specify the location of input data, the name of the output folder that will be created and the path of the sample sheet with *data_input_folder*, *output_folder* and *sample_sheet*  respectively.

* **Phenotypes**: here it is possible to specify the markers of interest into *pheno_list* (if not specified the complete list of phenotypes will be considered). This way only these markes will be considered for the analysis.

* **Descriptive**: here is possibile to specify parameters for descriptive section, as *raw* and/or *normalized* count evaluation and, when possible, statistical confrontation.

* **Map_plot**: here is possible to specify the markers to plot by *pheno_list* (if not specified the complete list of phenotypes will be considered).

* **Distance**: here are specified the parameters for distance calculation as *save_csv* if you want to save the distance values,  *pheno_list* (if not specified the complete list of phenotypes will be considered) if you want to calculate distance between specific markers, *pheno_from* and *pheno_to* to set specific couples of markers and distance direction, *plot_distance* to save distance curves and *save_csv_zscores* to save z-scores distance.

* **Image_match**: here is possibile to specify a sublist of markers to print on composite images (if not specified the complete list of phenotypes will be considered) via *pheno_list*. It is also possible to create interactive images, in addition to the static ones, by setting *interactive* option as true. Interactive graphs' layout can also be customize with the options *layout_marker_edge_col* and *layout_marker_size* to change respectively the edge color and the size of dots plotted on image and *layout_xsize* *layout_ysize* to customize the size of the image.

* **Distance_match**: here is possibile to specify a sublist of markers whose distances are to be printed on composite images via *pheno_list* (if not specified the complete list of phenotypes will be considered).

* **statistical_distance**: here is possibile to specify the markers for which you can perform the distance statsical analysis in *pheno_from* and *pheno_to*

* **Cluster**: here is possibile to specify the parameters for clustering analysis as *pheno_list* if you want to specify a restricted list of markers, *algo_method* to select the method to find optimal k value, *k* if you want to specify the maximum number of clusters, *cluster_method* to choose the clustering method (spectral, kmeans and k-prototype)

* **Cross_pcf**: here is possibile to specify the parameters for cross-PCF analysis as *radiiusOfInterest*, *anulusStep* and *anulusWidth* (for more info about this parameters check [Cross-PCF](https://www.cambridge.org/core/journals/biological-imaging/article/extended-correlation-functions-for-spatial-analysis-of-multiplex-imaging-data/FB677F0E100658E36725C5B4A3944EB7)). It is also possible to to specify a restricted list of markers through *pheno_list*, to plot a single image with all pcfs with *all_pcf* or to plot TCM maps on ROIs settign *on_roi*

* **Stats**: here is possible to set parameters for the statistical analysis like *sample_type* to define if the group is paired or unpaired and *p_adj* to set a p-value correction method.

#INSERIRE WORFLOW

### Input Folder Structure

The input folder for ALOA tool must be organized as follows:

```
input
├── img_match
│   ├── sbj001_[14146,53503]_composite_image.jpg
│   ├── ... 
│   └── sbjN_[13394,50883]_composite_image.jpg
│
├── raw_data
│   ├── sbj001
│   │   ├── sbj001_[14146,53503]_cell_seg_data.txt
│   │   ├── ...
│   │   └── sbj001_[17241,54367]_cell_seg_data.txt
│   ├── ...
│   └── sbjN
│       ├── sbjN_[13394,50883]_cell_seg_data.txt
│       ├── ...
│       └── sbjNT_[17130,56449]_cell_seg_data.txt
│
├── cellType_dict.tsv
└── sample_sheet.tsv

```
Where input is the main folder containing:
* raw_data: this folder contains as many subfolders as the number of subjects analyzed. Each subject folder contains the segmentation data for each one of the ROIs. <br>⚠️ This folder is required when -I and -D options are selected!
* img_match: this folder contains the ROIs on which plot markers and distances
* sample_sheet.tsv: this file must be compiled by the user with the subjects IDs and the belonging group. It's structure is shown below:

    | sbj_ID | Group|
    |----------------|----------------|
    |sbj001| Group1|
    |sbj002| Group1|
    |...| ...|
    |sbjN| GroupX|
* cellType.tsv: this file must be compiled by the user with the cell type name corresponding to each one of the markers of interest. <br>⚠️ This file is required when -p option is selected!

    | Phenotype | Cell_Type|
    |----------------|----------------|
    |FoxP3+| Regulatory T cell|
    |...| ...|
    |CD68+| Cytotoxic T cell|

### Output Folder Structure

All ALOA results are saved into a folder output whose path is specified in *config.json* file. The output folder is organized as follow:
```
output
├── Merged
├── Merged_clean
├── Maps_plot
├── Descriptive
├── Distance
├── Distance_Statistical
├── Clustering
├── Cross_PCF
├── Distance_match
├── Img_match
└── Log
```
All log files are saved into Log folder. The other output subfolders are named as the corresponding analysis and will be explained below.

### 1. Merge Cell Seg Data
In this section the cell_seg_data.txt files for each patient are merged into a single file through the integration of phenoptr [merge](https://akoyabio.github.io/phenoptrReports/articles/consolidation.html) function.

For the merge, two types of files are genetared:
- files where **negative cellulas** for markers of interest haven't been deleted (saved into *Merged* folder)
- files where **negative cellulas** for  markers of interest have been deleted (saved into *Merged_clean* folder)

Into *Merged* and *Merged_clean* folders as many folders as groups will be created. These subfolders will contain merged files for each patient of each group. Each merged file will be named as *patientnames.txt*

```
Merged
├── Group1
│   ├── Merge_cell_seg_data_sbj001.txt
│   ├── ...
│   └── Merge_cell_seg_data_sbjN.txt
├── ...
└── GroupN
    └── ...
Merged_clean
├── Group1
│   ├── Merge_cell_seg_data_sbj001.txt
│   ├── ...
│   └── Merge_cell_seg_data_sbjN.txt
├── ...
└── GroupN
    └── ...
```

Example:
```
python aloa.py -m
```

### 2. Map Plots
In this section a plot of all the markers of interest will be generate for each patient of each group as pdf images.

```
Maps_plot
├── Group1
│   ├── sbj001_Pheno_CD68+CD8+FoxP3+CK+.pdf
│   ├── ...
│   └── sbjN_Pheno_CD68+CD8+FoxP3+CK+.pdf
└── GroupN
    └── ...
```
Example:
```
python aloa.py -m -M
```

### 3. Descriptive data + statistical analysis
From merged file, this section produces raw and/or normalized markers counts, for each patient and for each group, saved into csv files (*csv* folder) and visualize through **barplot** figures saved into *Barplot* folder.

⚠️ Two different formulas are used for normalized counts in barplots (*Norm_count_patientname.csv*) and in boxplots (*all_norm_count_patientsname*). For more details see [functions ](./functions.md) section.


A *Descriptive* folder will be created to save all the results of this section.

<p align="center"><img src="Image_readme/subfolder.png" width=320></p>

From raw/normalized counts, if there are two or more groups, a statistical comparison is made to seek for significance difference between markers count. The comparison is shown through a **box plots** figure, with statistical annotation make though [TAP ](https://github.com/Discovery-Circle/tap)library, saved into *Box_Plot* folder

```
output_folder
Descriptive
├── Box Plots
│   ├── box_plot_comparison_Normalized.jpeg
│   └── box_plot_comparison_Raw.jpeg
├── Group1
│   ├── Bar_plot
│   │   ├── Bar_Plot_Normalized.jpeg
│   │   └── Bar_Plot_Raw.jpeg
│   └── csv
│       ├── Norm_count_sbj001.csv
│       ├── ...
│       ├── Raw_count_sbj001.csv
│       └── ...
├── ...
└── GroupN
    └── ...
            
```
Example:

```
python aloa.py -m -o
```
### 4. Distances Calculation

From merged file distances between markers are evaluated through phenoptr [distance](https://akoyabio.github.io/phenoptr/articles/computing_distances.html) function.

```
Distance
├── Group1
│   ├── sbj001_Distance.txt
│   ├── ...
│   └── sbjN_Distance.txt
├── ...
└── GroupN
    └── ...
```
Example:

```
python aloa.py -m -d
```
A statistical analysis can be also performed if two or more groups are reported.
A *Distance_Statistical* folder is created to save all the results of this section.

 <p align="center"><img src="Image_readme/distance_statistical.png" width=320> </p>

⚠️ Distance values are calculated as Z-score. For more details see For more details see [functions](./functions.md)


The results are: 
- **boxplot** figure,with statistical annotation make though [TAP ](https://github.com/Discovery-Circle/tap)library, saved in box_plot folder
- **csv** file with distance value for all patients for each group for a specific distance
- **distance_curve** figure (optional), saved into distance_curve folder, where are plotted the distributions of the distances of comparative groups for a specific couple of marker

```
Distance_Statistical
    ├── box_plot
    │   ├── distances_box_plot_CD68+_to_CD8+.png
    │   │   ...
    │   └── distances_box_plot_FoxP3+_to_CK+.png
    ├── csv
    │   ├── df_statistical_distance_CD68+_to_CD8+.csv
    │   │    ...
    │   └── df_statistical_distance_FoxP3+_to_CK+.csv
    ├── distance_curve
    │   ├── plot_statistical_distance_CD68+_to_CD8+.png
    │   │    ...
    │   └── plot_statistical_distance_FoxP3+_to_CK+.png
    └── summary_statistical.csv

```
Example:
```
python aloa.py -m -d -s
```

### 5. Clustering
From merged file clustering evaluation can be done. The output folder will be organized for each group as:
- **k optimal** folder containing all of the method curves images
- **cluster algorithm** folder containing *percentage* subfolder with all of the marker percentage in each cluster, *scatter_plot* containing a scatter for each subject and *stacked_barplot* with a stacked plot for each of the subject analyzed.
```
Clustering
├── Group1
│   ├── elbow_scores
│   │   │   ├── sbj001_[14146,53503].tiff
│   │   │   ├── ...
│   │   │   └── sbjN_[13394,50883].tiff
│   ├── kmeans
│   │   ├── percentage
│   │   │   ├── cluster_percentage_sbj001_[14146,53503].csv
│   │   │   ├── ...
│   │   │   └── cluster_percentage_sbjN_[13394,50883].csv
│   │   ├── scatter_plot
│   │   │   ├── sbj001_[14146,53503].tiff
│   │   │   ├── ...
│   │   │   └── sbjN__[13394,50883].tiff
│   │   └── stacked_barplot
│   │       ├── sbj001_[14146,53503].tiff
│   │       ├── ...
│   │       └── sbjN__[13394,50883].tiff
│   ├── prototype
│   │   ├── percentage
│   │   │   └── ...
│   │   ├── scatter_plot
│   │   │   └── ...
│   │   └── stacked_barplot
│   │       └── ...
│   └── spectral
│       ├── percentage
│       │   └── ...
│       ├── scatter_plot
│       │   └── ...
│       └── stacked_barplot
│           └── ...
├── ...
└── GroupN
    └── ...
```
Example:
```
python aloa.py -m -c
```

### 6. Cross-PCF
Cross-PCF is implemented by the integration of J. Bull cross-PCF [scripts](https://github.com/JABull1066/ExtendedCorrelationFunctions). For each group, a subfolder named as the subject ID is created and structured as:
- **Celltype couple** folder containing for each ROI and for each radius chosen three different images: cross-pcf curve, topographical map and point cloud plot.
- **summary** folder containing the csv files with all the cross-pcf values for each couple of celltypes.

If two or more groups are present, also a statistical analysis between cross-PCFs are evaluated and saved in **stats** folder for each radious.

```
Cross_PCF
├── Group1
│   ├── sbj001
│   │   ├── Macrophages-Cytotoxic_T_cell
│   │   │   ├── ROI_14146,53503
│   │   │   │   └── r_50
│   │   │   │       ├── PCF_function.tif
│   │   │   │       ├── TCM.tif
│   │   │   │       └── point_cloud.tif
│   │   │   ├── ...
│   │   │   └── Regulatory_T_cell-Cytotoxic_T_cell
│   │   │       └── ...
│   ├── ...
│   │   └── ...
│   ├── ...
│   └── summary
│       ├── Macrophages-Cytotoxic_T_cell.tsv
│       ├── ...
│       └── Regulatory_T_cell-Macrophages.tsv
│
├── ...
├── GroupN
│   └── ...
└── stats
    └── stat_analysis_r_50.tsv
```
Example:
```
python aloa.py -p
```
### 7. Imaging
The markers' position and their distances can be plot on single chosen ROIs images.
In particular, for marker positions the output folder will be structured as:
- **tif** images with markers plot on ROIs
- **Interactive_plots** folder containing the same images mentioned above but with the possibility to selecting and showing only specific markers by clicking on the correspondive legends dots.
```
Img_match
├── Interactive_plots
│   ├── sbj001_[14146,53503]_composite_image_match_CD68+CD8+FoxP3+CK+.html
│   │   ...
│   └── sbjN_[13394,50883]_composite_image_match_CD68+CD8+FoxP3+CK+.html
│
├── sbj001_[14146,53503]_composite_image_match_CD68+CD8+FoxP3+CK+.tif
│   ...
└── sbjN_[13394,50883]_composite_image_match_CD68+CD8+FoxP3+CK+.tif
```
Example:
```
python aloa.py -I
```
For distances plot a tiff image of the distances on ROI will be create for each couple of markers.
```
Distance_match
├── sbj001_[14146,53503]_composite_image_dist_match_CD68+CD8+_Nearest_CD68+_to_each_CD8+.tif
├── ...
└── sbjN_[13394,50883]_composite_image_dist_match_FoxP3+CK+_Nearest_FoxP3+_to_each_CK+.tif
```
Example:
```
python aloa.py -D
```