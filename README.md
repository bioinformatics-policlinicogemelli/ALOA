# ALOA (Analysis spatiaL prOfiling imAging)

## Introduction

ALOA is a useful bioinformatics tool for analyzing spatial data and imaging derived from Akoya PhenoimagerTM digital pathology workflow. For more informations [Unlock the Power of Spatial Biology with phenoptrReports](https://www.akoyabio.com/wp-content/uploads/2022/01/Spatial-Biology-with-phentoprReports-TechNote.pdf)

## Features
- Spatial **Data** Analysis 📈

    The first step is the merge of cell seg data files of each roi. This files are merged into a single file, for each patient. <br>With the option "", it's possibile to remove cells that don't have any of the phenotypes of interest (OTHER). <br>

    The results can be saved into:<br>
    - <em>output/Merged_clean</em><br>
                or
    - <em>output/Merged</em><br>



    From Merged or Merged clean file you can utilize 3 sections:
    - **Map Plot**, this section produces images to visualize spatial distribution of phenotypes of interest
    
 
    - **Description**, this section procudes **bar plots**, through which visualize phenotype count for each patient, and **box plots** through which compare phenotype count from groups. <br> For the count there is the possibility to calculate: <br>
            - **RAW** count         <br>
            - **NORMALIZED** count (inserire foto formule) <br>
    - **Distance**, this section calculated the distance from specific phenotypes or from the combination of all the phenotypes of the panel, starting from [phenoptr - Computing inter-cellular distances ](https://akoyabio.github.io/phenoptr/articles/computing_distances.html).<br> From distance data, it's possibile to make statistical analysis from groups, through box plot (default) and distribution curve plot (optional)

- Spatial **Imaging** Analysis <img src="images.jpeg" width=20 height=12>


## Installation


## Usage

