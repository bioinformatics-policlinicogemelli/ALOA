## In-depth functions
- **function_descriptive_analysis.py**

    The function_descriptive_analysis function is a Python script that allows you to make a phenotype description from the Merged_clean or Merged file and a statistical comparson if there are more than one group.<br>
    Through this function a parent folder is created, **descriptive**, containing just as many subfolders as there are groups present.<br>
    **Example**:<br> if we have 2 groups ("Group 2" and "Group 3")

    <img src="subfolder.png">


    Through this function you can:

    - Calculated the raw count of the phenotypes for each patient for each group,saved as *"Raw_count_patientID.csv"* into **csv** folder

    - Calculated the normalized count of the phenotypes for each patient for each group, saved as *"Norm_count_patientID.csv"* into **csv** folder<br>

    **++INSERIRE FORMULA NORMALIZZAZIONE++**<br>

    Starting from this count, a barplot is created to visualize the distribution of the phenotypes for each patient into a single group, sved the image into **Bar_plot** folder.

    **++INSERIRE IMMAGINE DI ESEMPIO DI UN BARPLOT++**

    If there are **more than one group**, this function makes a comparison between the group for the phenotypes of interest. Box plot are created for raw and/or normalized count and saved into Box Plot folder.<br>
    <br>
    <img src="boxplot_folder.png" height=300>
    <br>
    Statistical test is:
    - **Mann-Whitney** test if there are 2 groups
    - **Kruskal** test if there are more than 2 groups


    


   