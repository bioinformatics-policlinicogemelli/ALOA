import csv
import os
import re
import json
import pandas as pd
from scipy.stats import mannwhitneyu
import fitter
from loguru import logger
import pathlib

def split_and_calculate_means(filepath, pz_grado2, pz_grado3):
    # Estrazione dei nomi C_1 e C_2 dal nome del file
    match = re.match(r'.*/risultati_(.*?)_(.*?)\.csv', filepath)
    if match:
        C_1 = match.group(1)
        C_2 = match.group(2)

        # Lettura del file CSV di input e divisione delle righe nei due gruppi
        data = pd.read_csv(filepath)

        # Divisione in gruppi grado2 e grado3
        grado2 = data[data['Paziente'].isin(pz_grado2)]
        grado3 = data[data['Paziente'].isin(pz_grado3)]

        # Calcolo delle medie
        media_C1_grado2 = grado2['Conteggio_C1'].mean()
        media_C2_grado2 = grado2['Conteggio_C2'].mean()
        media_C1_grado3 = grado3['Conteggio_C1'].mean()
        media_C2_grado3 = grado3['Conteggio_C2'].mean()
        PCF_r_grado2 = grado2['PCF_r'].mean()
        PCF_r_grado3 = grado3['PCF_r'].mean()

        # Test di Mann-Whitney
        mw_stat, mw_p_value = mannwhitneyu(grado2['PCF_r'], grado3['PCF_r'])
        #print(stat, p_value)

        # Fit della distribuzione per PCF_r_grado2 e PCF_r_grado3
        best_fit_grado2 = fit_distribution(grado2['PCF_r'])
        best_fit_grado3 = fit_distribution(grado3['PCF_r'])

        print(f"Best fit distribution for Grado 2: {best_fit_grado2}")
        print(f"Best fit distribution for Grado 3: {best_fit_grado3}")


        return (C_1, C_2, media_C1_grado2, media_C2_grado2, media_C1_grado3, media_C2_grado3, PCF_r_grado2, PCF_r_grado3, mw_stat, mw_p_value, best_fit_grado2, best_fit_grado3)

def fit_distribution(data):
    f = fitter.Fitter(data, distributions=['norm', 'expon', 'uniform', 'gamma', 'beta', 'lognorm'])
    f.fit()
    return f.get_best(method='sumsquare_error')



def main():
    
    print("\n################################# PCF ANALYSIS #################################\n")
    
    logger.info("Start pcf analysis: This step will perform a pcf analysis on single ROIs\n")

    logger.info("Reading configuration file")
    with open("config.json") as f:
        data=json.load(f)
    
    input_folder=os.path.join(data["Paths"]["data_input_folder"],"img_match")
   
    output_path = os.path.join(data["Paths"]["output_folder"],"PCF")
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)

    # Definizione dei codici per i due gruppi
    pz_grado2 = ["10647", "15397", "4315", "5808"]
    pz_grado3 = ["07B07364", "13B12286", "14B03092", "14B26993", "8033K1", "9998"]
    
    with open("/home/beatrice/Documenti/università/MAGISTRALE/tirocinio Gemelli/ExtendedCorrelationFunctions/script Bea/config.json") as f:
        data = json.load(f)
    
    results = []
    directory = data["analisi_PCF"]["path"]
    output_file = "medie_result_no_14b04959.csv" 

    # Itera su tutti i file nella directory
    for filename in os.listdir(directory):
        if filename.startswith('risultati') and filename.endswith('.csv'):
            filepath = os.path.join(directory, filename)
            result = split_and_calculate_means(filepath, pz_grado2, pz_grado3)
            results.append(result)
    try:
        print(f"Writing results to: {output_file}")
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['C_1', 'C_2', 'media_C1_grado2', 'media_C2_grado2', 'media_C1_grado3', 'media_C2_grado3', 'PCF_r_grado2', 'PCF_r_grado3', 'stat_MW', 'p_value_MW', 'best_fit_grado2', "best_fit_grado3"])
            writer.writerows(results)
        print("Calcolo delle medie e test di Mann-Whitney completato. I risultati sono stati salvati in 'medie_result.csv'.")
    
    except Exception as e:
        print(f"Error writing results: {e}")

    # # Scrivi i risultati su un nuovo file CSV
    # with open('medie_result.csv', 'w', newline='') as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerow(['C_1', 'C_2', 'media_C1_grado2', 'media_C2_grado2', 'media_C1_grado3', 'media_C2_grado3', 'PCF_r_grado2', 'PCF_r_grado3', 'stat', 'p_value'])
    #     writer.writerows(results)


if __name__ == "__main__":
    main()
