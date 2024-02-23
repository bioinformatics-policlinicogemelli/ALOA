# ALOA (Analysis spatiaL prOfiling imAging)

## Introduction

ALOA is a useful bioinformatics tool for analyzing spatial data and imaging derived from Akoya PhenoimagerTM digital pathology workflow ([Unlock the Power of Spatial Biology with phenoptrReports](https://www.akoyabio.com/wp-content/uploads/2022/01/Spatial-Biology-with-phentoprReports-TechNote.pdf))

## Features
- Spatial **Data** Analysis 📈

    The first step is the merge of cell seg data files of each roi. This files are merged into a single file, for each patient. <br>With the option "", it's possibile to remove cells that don't have any of the phenotypes of interest (OTHER). <br>

    The results can be saved into:<br>
    - <em>output/Merged_clean</em><br>
                or
    - <em>output/Merged</em><br>



    From *Merged* or *Merged_clean* files you can utilize 3 sections:
    - **Map Plot**, this section produces images to visualize spatial distribution of markers of interest
    
 
    - **Description**, this section procudes **bar plots**, through which visualize markers count for each patient, and **box plots** through which compare markers count from groups. <br> For the count there is the possibility to calculate: <br>
            - **RAW** count         <br>
            - **NORMALIZED** count <br>
    - **Distance**, this section calculated the distance from specific markers or from the combination of all the markers of the panel, starting from [phenoptr - Computing inter-cellular distances ](https://akoyabio.github.io/phenoptr/articles/computing_distances.html).<br> From distance data, it's possibile to make statistical analysis from groups, through box plot (default) and distribution curve plot (optional)

- Spatial **Imaging** Analysis <img src="Image_redme/images.jpeg" width=20 height=12>


## Installation


## Usage

#### Launch ALOA main

These are the options that can be set for this block:

| Options | Input | Type | Required
|----------------|----------------| :---:| :---:|
|-m <br> --merge| <p align="justify">merge datasets from ROIs of the same patient| boolean | No
|-c <br> --clean| <p align="justify">remove the all Phenotype-negative cells| boolean | No
|-M <br> --maps| <p align="justify">create maps plot| boolean | No
|-d <br> --distance| <p align="justify">distance evaluation between phenotypes| boolean | No
|-s <br> --stats| <p align="justify">create distance stats| boolean | No
|-o <br> --overview| <p align="justify">create data overview | boolean | No
|-C <br> --clustering| <p align="justify">cluster data| boolean | No
|-I <br> --clustering| <p align="justify">create phenotypes image match| boolean | No
|-D <br> --clustering| <p align="justify">create phenotypes distance match| boolean | No
|-w <br> --clustering| <p align="justify">create data overview| boolean | No
