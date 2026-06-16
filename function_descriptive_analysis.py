#Copyright 2024 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#    http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

from loguru import logger
import pathlib
import os
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math
import tap
import glob

# --- New imports used ONLY for the additional statistics table/annotation.
# The descriptive counts, normalization and original tap.plot_stats calls are unchanged.
# scipy serve solo per calcolare esplicitamente i p-value da salvare nel CSV.
# Se non è installato, il resto dello script continua comunque a funzionare.
try:
    from scipy import stats
except Exception:
    stats = None

# statsmodels serve solo per calcolare il padj con metodi standard
# come bonferroni, holm, fdr_bh, fdr_by.
# Anche questo è opzionale: sotto ci sono dei fallback per i metodi più comuni.
try:
    from statsmodels.stats.multitest import multipletests
except Exception:
    multipletests = None


#******************************************************
# Conta, per ogni gruppo e per ogni paziente, quante cellule appartengono
# ai fenotipi indicati nel config in Phenotypes -> pheno_list.
# Questa funzione è invariata rispetto alla logica originale.
def raw_count_cells(PATH_MERGE_FOLDER,list_pheno):
    logger.info("Start raw count calculation:")
    dict_global_info = {}
    for _d in [f for f in os.listdir(PATH_MERGE_FOLDER) if not f.startswith('.')]:
        sub_directory=os.path.join(PATH_MERGE_FOLDER,_d)
        if not os.listdir(sub_directory):
            logger.error(f"Directory {_d} is empty")
            exit()
        for f in [f for f in os.listdir(sub_directory) if not f.startswith('.')]:
            _id_paz=f.replace("Merge_cell_seg_data_clean_","").replace(".txt","")
            logger.info(f"Subject: {_id_paz}")
            dict_global_info.setdefault(_d,{})
            dict_global_info[_d].setdefault(_id_paz,{})
            _file=os.path.join(sub_directory,f)
            _data=pd.read_csv(_file, sep="\t")
            _data["Pheno"]=list(map(lambda x: x.replace("+","+,"),_data["Pheno"]))
            _data["Pheno"]=list(map(lambda x: x[:-1],_data["Pheno"]))
            logger.info(f"Reading File {_file}")
            _total_cells=len(_data)
            dict_global_info[_d][_id_paz]["Total_Cells"]=_total_cells
            _data_groupped = _data.groupby(["Pheno"])["Pheno"].count()
            _data_dict = {frozenset(k.replace("OTHER,", "").replace(",OTHER", "").split(",")): v for k, v in _data_groupped.to_dict().items()}
            for main_pheno in list_pheno:
                dict_global_info[_d][_id_paz][f"{main_pheno}"] = _data_dict.get(frozenset(main_pheno.split(",")), 0)
    logger.info("End raw count calculation")
    return dict_global_info 

#******************************************************
def calculate_mean_group_cells(dictionary_raw_count):
    count_mean_group_cells={}
    for _k,_v in dictionary_raw_count.items():
        _temp=[]
        count_mean_group_cells.setdefault(_k,0)
        for _,_v2 in _v.items():
            _temp.append(_v2["Total_Cells"])
        count_mean_group_cells[_k]=round(np.mean(_temp),2)
    return count_mean_group_cells

#******************************************************
def normalized_count_cells(dictionary_raw_count,dictionary_mean_count):
    logger.info("Start normalized count calculation")
    dict_normalized_count_cells={}
    for group, patients in dictionary_raw_count.items():
        dict_normalized_count_cells[group]={}
        mean_group=dictionary_mean_count[group]
        for patient, counters in patients.items():
            logger.info(f"Subject: {patient}")
            dict_normalized_count_cells[group][patient]={}
            total_cells=counters["Total_Cells"]
            dict_normalized_count_cells[group][patient]["Total_Cells"]=total_cells
            for pheno, value in counters.items():
                if pheno=="Total_Cells": continue
                _new_value=(value/total_cells)*mean_group
                dict_normalized_count_cells[group][patient][pheno]=_new_value
    logger.info("End normalized count calculation")
    return dict_normalized_count_cells

#******************************************************
def create_output_dir(path_output,groups_name):
    for _g in groups_name:
        temp_folder=os.path.join(path_output,_g)
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
            logger.info(f" Created folder {temp_folder}")

#******************************************************
def create_summury_file(path_output_results,dictionary_count,type_data):
    logger.info("writing csv file with count info for each patient")
    for group,patients in dictionary_count.items():
        logger.info(f"Group: {group}")
        dire = os.path.join(path_output_results,group,"csv")
        pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        for patient,phenos in patients.items():
            logger.info(f"Subject: {patient}")
            with open(os.path.join(dire, f"{type_data}_count_{patient}.csv"),"w") as f:
                f.write(f"Patient\tPheno\tCount_{type_data}\n")
                logger.info(f"Start Writing csv file...")
                for p in phenos:
                    f.write(f'{patient}\t{p}\t{phenos[p]}\n')
                logger.info(f"End Writing {type_data}_count_{patient}.csv for {patient}")

#******************************************************
def create_norm_all_file(path_output_results,dictionary_norm_count):
    logger.info("Start normalized count calculatin on all groups")
    for group,patients in dictionary_norm_count.items():
        dire = os.path.join(path_output_results,group,"csv")
        pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        for patient,phenos in patients.items():
            with open(os.path.join(dire, "all_norm_count_"+patient+".csv"),'w') as f:
                f.write(f"Patient\tPheno\tCount_all_Norm\n")
                logger.info(f"Writing normalized count file on all groups for patient {patient}")
                for p in phenos:
                    f.write(f'{patient}\t{p}\t{phenos[p]}\n')
    logger.info("End normalized count calculation on all groups")

#******************************************************
# Unisce i CSV per-paziente in una tabella unica wide:
# una riga = un paziente, colonne = fenotipi.
# Utile per avere merged_raw.csv, merged_norm.csv, merged_all_norm.csv.
def merge_counts(path_output_results, pattern, output_filename, phenos_list):
    all_files = []
    for group_dir in [d for d in os.listdir(path_output_results) if not d.startswith('.')]:
        csv_dir = os.path.join(path_output_results, group_dir, "csv")
        files = glob.glob(os.path.join(csv_dir, pattern))
        all_files.extend(files)

    if not all_files:
        logger.warning(f"No files found for pattern {pattern}")
        return

    # Load all files
# Load all files and reconstruct Group from folder structure
    df_list = []
    for f in all_files:
        group = pathlib.Path(f).parents[2].name
        tmp = pd.read_csv(f, sep="\t")
        tmp["Group"] = group
        df_list.append(tmp)

    df = pd.concat(df_list, ignore_index=True)


    # Infer the count column automatically
    count_col = [c for c in df.columns if c.startswith("Count_") or c.startswith("Count_all")][0]

    # Pivot to WIDE format
    df_wide = df.pivot_table(
        index=["Group", "Patient"],
        columns="Pheno",
        values=count_col,
        aggfunc="first"
    ).reset_index()

    # Ensure all phenotypes exist as columns (even if missing)
    for ph in phenos_list:
        if ph not in df_wide.columns:
            df_wide[ph] = 0

    if "Total_Cells" not in df_wide.columns:
        df_wide["Total_Cells"] = 0

    # Reorder columns
    ordered_cols = ["Group", "Patient"] + phenos_list + ["Total_Cells"]
    existing = [c for c in ordered_cols if c in df_wide.columns]
    df_wide = df_wide[existing]

    # Save final wide file
    merged_path = os.path.join(path_output_results, output_filename)
    df_wide.to_csv(merged_path, sep="\t", index=False)

    logger.info(f"✅ WIDE merged file saved as {merged_path}")

#******************************************************
# Crea i CSV singoli e poi li unisce nei file merged_*.csv.
# Questa funzione mantiene la struttura di output originale.
def create_and_merge_counts(path_output_results, dict_raw_count, phenos_list,
                            dict_norm_count=None, dict_norm_all_count=None):

    create_summury_file(path_output_results, dict_raw_count, "Raw")
    merge_counts(path_output_results, "Raw_count_*.csv", "merged_raw.csv", phenos_list)

    if dict_norm_count is not None:
        create_summury_file(path_output_results, dict_norm_count, "Norm")
        merge_counts(path_output_results, "Norm_count_*.csv", "merged_norm.csv", phenos_list)

    if dict_norm_all_count is not None:
        create_norm_all_file(path_output_results, dict_norm_all_count)
        merge_counts(path_output_results, "all_norm_count_*.csv", "merged_all_norm.csv", phenos_list)

#******************************************************
# Crea i barplot dei conteggi per gruppo/paziente.
# Non è stata modificata la logica grafica originale.
def bar_plot(path_output_result,dict_data,type_data):
    logger.info(f"Barplot creation for {type_data} count")
    for group,data in dict_data.items():
        dire = os.path.join(path_output_result,group,"Bar_plot")
        pathlib.Path(dire).mkdir(parents=True, exist_ok=True)
        logger.info(f"Group: {group}")
        max_val_y=0
        fig=go.Figure()
        for patient,targets in data.items():
            logger.info(f"Subject: {patient}")
            sorted_target=dict(sorted(targets.items(),key=lambda x:x[0]))
            target=list(sorted_target.keys())
            target.remove("Total_Cells")
            values=list(sorted_target.values())[:-1]
            max_val=max(values)
            if max_val>max_val_y:
                max_val_y=max_val
            fig.add_trace(go.Bar(x=target,y=values, name=patient))
        fig.update_layout(barmode='group',title_text=f"Phenotypes {type_data}-{group}",xaxis=dict(
            title="Phenotypes",showticklabels=True),yaxis=dict(title=f'log({type_data} Count)',showticklabels=True,type="log"),legend_title_text="Patients")
        fig.update_yaxes(range=[0,(math.log(max_val_y,10)+1.0)])
        fig.write_image(os.path.join(dire,"Bar_Plot_"+type_data+".jpeg"),scale=6)

#******************************************************
def calculate_mean_total_groups(dictionary_raw_count):
    _count_total=[]
    for _k,_v in dictionary_raw_count.items():
        for _k2,_v2 in _v.items():
            _count_total.append(_v2["Total_Cells"])
    mean_cells=np.mean(_count_total)
    return mean_cells

#******************************************************
# Normalizzazione rispetto alla media del numero totale di cellule su tutti i gruppi.
# Questa è la normalizzazione usata poi per i boxplot "Normalized".
def normalized_count_on_all_groups(dictionary_raw_count,mean_all_groups):
    logger.info("Starting normalized count on all groups")
    dict_normalized_count_cells={}
    for group, patients in dictionary_raw_count.items():
        logger.info(f"Group: {group}")
        dict_normalized_count_cells[group]={}
        for patient, counters in patients.items():
            logger.info(f"Subject: {patient}")
            dict_normalized_count_cells[group][patient]={}
            total_cells=counters["Total_Cells"]
            dict_normalized_count_cells[group][patient]["Total_Cells"]=total_cells
            for pheno, value in counters.items():
                if pheno == "Total_Cells": continue
                _new_value=(value/total_cells)*mean_all_groups
                dict_normalized_count_cells[group][patient][pheno]=_new_value
    logger.info("End normalized count on all groups")
    return dict_normalized_count_cells

#******************************************************
# Converte il dizionario dei conteggi in un DataFrame lungo
# con colonne group, pheno, Count, cioè il formato richiesto dai boxplot.
def prepare_data_box_plot(dictionary_count):
    _data_all=pd.DataFrame()
    _data_norm=[]
    for _k,_v in dictionary_count.items():
        for _,_v2 in _v.items():
            for _t in _v2.keys():
                if _t=="Total_Cells": continue
                _data_norm.append([_k,_t,_v2[_t]])
    _data_all=pd.DataFrame(_data_norm,columns=["group","pheno","Count"])
    return _data_all

#******************************************************
def create_output_box_plot_dir(path_output):
    temp_folder=os.path.join(path_output,"Box Plots")
    pathlib.Path(temp_folder).mkdir(parents=True, exist_ok=True)


#******************************************************
# ======================================================================
# FUNZIONI AGGIUNTE PER GESTIRE PVALUE/PADJ DA CONFIG
# ======================================================================
# Da qui inizia la parte aggiunta rispetto allo script originale.
# Serve a: leggere il metodo di correzione, calcolare pvalue/padj,
# salvare un CSV statistico e, se richiesto, annotare i boxplot.
def _normalize_p_adjust(p_adjust):
    """
    Normalize the p-value correction name coming from the config.

    Accepted values in config:
      - "bonferroni" or "bonf"
      - "holm"
      - "fdr_bh", "fdr", "bh", "benjamini-hochberg"
      - "fdr_by", "by", "benjamini-yekutieli"
      - "none", "", None  -> no correction; padj = pvalue

    This function is used for the CSV/statistical annotation. The same normalized
    value is also passed to tap.plot_stats through type_correction, exactly where
    the original script already passed p_adjust.
    """
    if p_adjust is None:
        return None

    p_adjust = str(p_adjust).strip().lower()
    if p_adjust in ["", "none", "no", "false", "null", "nan"]:
        return None
    if p_adjust in ["bonf", "bonferroni"]:
        return "bonferroni"
    if p_adjust in ["holm", "holm-bonferroni"]:
        return "holm"
    if p_adjust in ["fdr", "bh", "fdr_bh", "benjamini-hochberg"]:
        return "fdr_bh"
    if p_adjust in ["by", "fdr_by", "benjamini-yekutieli"]:
        return "fdr_by"

    logger.warning(f"Unknown p_adj='{p_adjust}'. It will be passed as-is to tap.plot_stats, but CSV padj will equal pvalue.")
    return p_adjust

#******************************************************
def _adjust_pvalues(pvalues, p_adjust):
    """
    Calculate adjusted p-values for the CSV and for the text annotation.

    Important: in the original script, the correction for the plot was delegated
    to tap.plot_stats using type_correction=p_adjust. Here we add an explicit
    padj column so that the output is inspectable and filterable.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    p_adj = np.full_like(pvalues, np.nan, dtype=float)

    valid = ~np.isnan(pvalues)
    if valid.sum() == 0:
        return p_adj

    method = _normalize_p_adjust(p_adjust)

    # No correction: padj is equal to the raw p-value.
    if method is None:
        p_adj[valid] = pvalues[valid]
        return p_adj

    # Preferred implementation, if statsmodels is available.
    if multipletests is not None:
        try:
            _, corrected, _, _ = multipletests(pvalues[valid], method=method)
            p_adj[valid] = corrected
            return p_adj
        except Exception as e:
            logger.warning(f"Could not adjust p-values with method '{method}' using statsmodels: {e}")

    # Fallback implementations for the most common methods.
    pv = pvalues[valid]
    m = len(pv)

    if method == "bonferroni":
        corrected = np.minimum(pv * m, 1.0)
    elif method == "holm":
        order = np.argsort(pv)
        sorted_p = pv[order]
        sorted_adj = np.maximum.accumulate((m - np.arange(m)) * sorted_p)
        sorted_adj = np.minimum(sorted_adj, 1.0)
        corrected = np.empty_like(sorted_adj)
        corrected[order] = sorted_adj
    elif method == "fdr_bh":
        order = np.argsort(pv)
        sorted_p = pv[order]
        sorted_adj = sorted_p * m / (np.arange(m) + 1)
        sorted_adj = np.minimum.accumulate(sorted_adj[::-1])[::-1]
        sorted_adj = np.minimum(sorted_adj, 1.0)
        corrected = np.empty_like(sorted_adj)
        corrected[order] = sorted_adj
    elif method == "fdr_by":
        c_m = np.sum(1 / np.arange(1, m + 1))
        order = np.argsort(pv)
        sorted_p = pv[order]
        sorted_adj = sorted_p * m * c_m / (np.arange(m) + 1)
        sorted_adj = np.minimum.accumulate(sorted_adj[::-1])[::-1]
        sorted_adj = np.minimum(sorted_adj, 1.0)
        corrected = np.empty_like(sorted_adj)
        corrected[order] = sorted_adj
    else:
        logger.warning(f"No fallback available for p_adj='{method}'. CSV padj will equal pvalue.")
        corrected = pv

    p_adj[valid] = corrected
    return p_adj

#******************************************************
def _get_stats_options(data):
    """
    Read the statistics options from the config.

    Original config already had:
      Stats -> sample_type
      Stats -> p_adj

    New optional config keys:
      Stats -> p_adj_threshold          e.g. 0.05 or 0.1
      Stats -> save_stats_csv           true/false
      Stats -> show_pvalues_on_boxplot  true/false
    """
    stats_config = data.get("Stats", {})

    # Metodo di correzione letto dal config, ad esempio:
    # bonferroni, holm, fdr_bh, fdr_by, none.
    p_adjust = _normalize_p_adjust(stats_config.get("p_adj", None))
    # Soglia di significatività: 0.05 per analisi standard,
    # 0.1 se vuoi una soglia più esplorativa.
    p_threshold = float(stats_config.get("p_adj_threshold", 0.05))

    # Se true, salva statistics_Raw.csv / statistics_Normalized.csv.
    save_stats_csv = bool(stats_config.get("save_stats_csv", True))

    # Se true, aggiunge una riga sotto ai boxplot con pvalue e padj.
    show_pvalues_on_boxplot = bool(stats_config.get("show_pvalues_on_boxplot", True))

    # sample_type = "paired" usa Wilcoxon; altrimenti Mann-Whitney per 2 gruppi.
    test = stats_config.get("sample_type", "")

    return p_adjust, p_threshold, save_stats_csv, show_pvalues_on_boxplot, test

#******************************************************
def _safe_format_pvalue(pvalue):
    """Format p-values in a compact and readable way for plot annotation."""
    if pd.isna(pvalue):
        return "NA"
    if pvalue < 1e-4:
        return f"{pvalue:.2e}"
    return f"{pvalue:.4f}"

#******************************************************
def _compute_stats_table(data_all, hue_order, labels, p_adjust, test, p_threshold):
    """
    Compute a transparent statistics table with raw p-value and adjusted p-value.

    This DOES NOT replace tap.plot_stats. It only creates an explicit CSV and the
    optional text that is appended under the image.

    For 2 groups:
      - Mann-Whitney if sample_type != "paired"
      - Wilcoxon if sample_type == "paired" and the two vectors have same length

    For >2 groups:
      - Kruskal-Wallis global test per phenotype
      - If pairwise post-hoc p-values are needed, tap.plot_stats still handles the
        Dunn plot as in the original script.
    """
    rows = []

    if stats is None:
        logger.warning("scipy is not available: statistics CSV will contain NA p-values.")

    # Il test viene calcolato separatamente per ogni fenotipo.
    for pheno in labels:
        df_pheno = data_all[data_all["pheno"] == pheno].copy()

        # Caso più comune: confronto tra due gruppi, ad esempio HER2 vs TNBC.
        if len(hue_order) == 2:
            g1, g2 = hue_order[0], hue_order[1]
            v1 = df_pheno[df_pheno["group"] == g1]["Count"].dropna().astype(float).values
            v2 = df_pheno[df_pheno["group"] == g2]["Count"].dropna().astype(float).values

            p_value = np.nan
            test_name = "wilcoxon" if test == "paired" else "mannwhitneyu"

            if stats is not None and len(v1) > 0 and len(v2) > 0:
                try:
                    if test == "paired":
                        if len(v1) == len(v2):
                            p_value = stats.wilcoxon(v1, v2).pvalue
                        else:
                            logger.warning(
                                f"Wilcoxon paired requested for {pheno}, but group sizes differ "
                                f"({g1}: {len(v1)}, {g2}: {len(v2)}). pvalue set to NA in CSV."
                            )
                    else:
                        p_value = stats.mannwhitneyu(v1, v2, alternative="two-sided").pvalue
                except Exception as e:
                    logger.warning(f"Could not calculate pvalue for {pheno} ({g1} vs {g2}): {e}")

            rows.append({
                "pheno": pheno,
                "test": test_name,
                "comparison": f"{g1}_vs_{g2}",
                "group_1": g1,
                "group_2": g2,
                "n_group_1": len(v1),
                "n_group_2": len(v2),
                "mean_group_1": np.mean(v1) if len(v1) else np.nan,
                "mean_group_2": np.mean(v2) if len(v2) else np.nan,
                "median_group_1": np.median(v1) if len(v1) else np.nan,
                "median_group_2": np.median(v2) if len(v2) else np.nan,
                "p_value": p_value,
            })

        # Se ci sono più di due gruppi viene salvato un test globale
        # Kruskal-Wallis per ogni fenotipo. Il post-hoc Dunn rimane gestito
        # da tap.plot_stats, come nello script originale.
        elif len(hue_order) > 2:
            values_by_group = []
            for g in hue_order:
                vals = df_pheno[df_pheno["group"] == g]["Count"].dropna().astype(float).values
                if len(vals) > 0:
                    values_by_group.append(vals)

            p_value = np.nan
            if stats is not None and len(values_by_group) >= 2:
                try:
                    p_value = stats.kruskal(*values_by_group).pvalue
                except Exception as e:
                    logger.warning(f"Could not calculate Kruskal-Wallis pvalue for {pheno}: {e}")

            rows.append({
                "pheno": pheno,
                "test": "kruskal_wallis_global",
                "comparison": "global",
                "group_1": "all_groups",
                "group_2": "",
                "n_group_1": int(df_pheno["Count"].notna().sum()),
                "n_group_2": np.nan,
                "mean_group_1": df_pheno["Count"].mean(),
                "mean_group_2": np.nan,
                "median_group_1": df_pheno["Count"].median(),
                "median_group_2": np.nan,
                "p_value": p_value,
            })

    df_stats = pd.DataFrame(rows)

    if len(df_stats) == 0:
        return df_stats

    # Qui viene calcolato il padj in modo esplicito.
    # Questa è la parte che prima non esisteva come colonna salvata.
    df_stats["p_adj"] = _adjust_pvalues(df_stats["p_value"].values, p_adjust)

    # Metodo e soglia vengono salvati nel CSV per tracciabilità.
    df_stats["p_adj_method"] = _normalize_p_adjust(p_adjust) if _normalize_p_adjust(p_adjust) is not None else "none"
    df_stats["p_adj_threshold"] = p_threshold

    # La significatività viene decisa solo in base al padj e alla soglia del config.
    df_stats["significant"] = df_stats["p_adj"] <= p_threshold

    return df_stats

#******************************************************
def _append_stats_to_image(filename, stats_table, max_lines=8):
    """
    Append pvalue/padj text below the JPEG generated by tap.plot_stats.

    This keeps the original plotting function unchanged and only adds a white
    strip at the bottom of the saved image. If PIL is unavailable or anything
    fails, the script keeps the original image and logs a warning.
    """
    if stats_table is None or len(stats_table) == 0:
        return

    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception as e:
        logger.warning(f"PIL not available; pvalue/padj annotation skipped for {filename}: {e}")
        return

    try:
        img = Image.open(filename).convert("RGB")
        draw_probe = ImageDraw.Draw(img)

        # Testo che viene aggiunto sotto all'immagine:
        # p = p-value grezzo; padj = p-value corretto.
        lines = []
        for _, row in stats_table.head(max_lines).iterrows():
            sig = "significant" if bool(row.get("significant", False)) else "not significant"
            lines.append(
                f"{row['pheno']} | {row['comparison']} | "
                f"p={_safe_format_pvalue(row['p_value'])}; "
                f"padj={_safe_format_pvalue(row['p_adj'])} "
                f"({row['p_adj_method']}, threshold={row['p_adj_threshold']}) | {sig}"
            )

        if len(stats_table) > max_lines:
            lines.append(f"... {len(stats_table) - max_lines} additional comparison(s) in the statistics CSV")

        font = ImageFont.load_default()
        line_height = 18
        padding = 12
        extra_h = padding * 2 + line_height * len(lines)

        # Crea una fascia bianca in basso e ci scrive le statistiche.
        # L'immagine del boxplot originale non viene ridisegnata.
        new_img = Image.new("RGB", (img.width, img.height + extra_h), "white")
        new_img.paste(img, (0, 0))
        draw = ImageDraw.Draw(new_img)

        y = img.height + padding
        for line in lines:
            draw.text((padding, y), line, fill="black", font=font)
            y += line_height

        new_img.save(filename)
        logger.info(f"Added pvalue/padj annotation to {filename}")

    except Exception as e:
        logger.warning(f"Could not append pvalue/padj annotation to {filename}: {e}")

#******************************************************
# Crea i boxplot comparativi.
# La funzione originale è stata estesa solo per ricevere:
# - soglia padj
# - salvataggio CSV statistiche
# - annotazione opzionale pvalue/padj sul boxplot
def create_comparison_box_plot(path_output_result, data_all, p_adjust, test, type_data,
                               p_threshold=0.05, save_stats_csv=True,
                               show_pvalues_on_boxplot=True):
    logger.info(f"Box Plot creation on {type_data} count")
    x="pheno"
    y="Count"
    hue="group"
    hue_order=list(data_all["group"].unique())
    labels=sorted(list(data_all["pheno"].unique()))
    base_dir=os.path.join(path_output_result,"Box Plots",type_data)
    pathlib.Path(base_dir).mkdir(parents=True, exist_ok=True)

    # New: explicit statistics table with pvalue, padj and significance threshold.
    # The original tap.plot_stats call below is still used for the actual boxplot.
    # Calcola la tabella statistica una volta sola per tutti i fenotipi.
    # Questa tabella viene poi salvata e usata per annotare i plot.
    stats_table = _compute_stats_table(data_all, hue_order, labels, p_adjust, test, p_threshold)
    if save_stats_csv and len(stats_table) > 0:
        stats_csv = os.path.join(base_dir, f"statistics_{type_data}.csv")
        stats_table.to_csv(stats_csv, sep="\t", index=False)
        logger.info(f"Statistics table saved as {stats_csv}")

    # Mega plot
    filename=os.path.join(base_dir,f"box_plot_comparison_{type_data}.jpeg")
    _run_stat_test(data_all, labels, hue_order, p_adjust, test, x, y, hue, filename, type_data,
                   stats_table=stats_table, show_pvalues_on_boxplot=show_pvalues_on_boxplot)
    # Single plots
    for pheno in labels:
        logger.info(f"Creating single boxplot for phenotype: {pheno}")
        df_single=data_all[data_all["pheno"]==pheno].copy()
        filename_single=os.path.join(base_dir,f"box_plot_{pheno}_{type_data}.jpeg")
        stats_single = stats_table[stats_table["pheno"] == pheno].copy() if len(stats_table) > 0 else None
        _run_stat_test(df_single, [pheno], hue_order, p_adjust, test, x, y, hue, filename_single, type_data,
                       stats_table=stats_single, show_pvalues_on_boxplot=show_pvalues_on_boxplot)

#******************************************************
# Wrapper della chiamata originale a tap.plot_stats.
# Il comportamento statistico/grafico principale resta quello di tap.plot_stats;
# in più, alla fine, può aggiungere il testo pvalue/padj sotto il JPEG.
def _run_stat_test(data, labels, hue_order, p_adjust, test, x, y, hue, filename, type_data,
                   stats_table=None, show_pvalues_on_boxplot=True):
    kwargs_common={"width":4000,"height":1000,"title":f"Comparison {type_data} Count","log_y":True,
                   "labels":{"pheno":"Phenotypes","value":f"log({type_data} Counts)","group":"Group"}}

    # Same role as the original p_adjust variable, but normalized from config.
    p_adjust_for_tap = _normalize_p_adjust(p_adjust)

    try:
        if len(hue_order)==2 and test!="paired":
            logger.info(f"Mann-Whitney Test → {filename}")
            tap.plot_stats(data, x, y, order=labels, filename=filename,
                           type_correction=p_adjust_for_tap, export_size=(1400,950,3),
                           subcategory=hue, kwargs=kwargs_common)
        elif len(hue_order)==2 and test=="paired":
            logger.info(f"Wilcoxon Test → {filename}")
            # Chiamata originale: per 2 gruppi paired usa Wilcoxon.
            tap.plot_stats(data, x, y, order=labels, filename=filename,
                           type_test="wilcoxon", type_correction=p_adjust_for_tap, export_size=(1400,950,3),
                           subcategory=hue, kwargs=kwargs_common)
        elif len(hue_order)>2:
            logger.info(f"Kruskal-Wallis Test → {filename}")
            if p_adjust_for_tap is None: p_adjust_for_tap="bonferroni"
            # Chiamata originale: con >2 gruppi usa Dunn/post-hoc tramite tap.plot_stats.
            tap.plot_stats(data, x, y, type_test="dunn", order=labels,
                           subcategory=hue, filename=filename,
                           type_correction=p_adjust_for_tap, export_size=(1400,950,3),
                           kwargs=kwargs_common)
        # New: append pvalue/padj text below the image, only if requested in config.
        if show_pvalues_on_boxplot:
            _append_stats_to_image(filename, stats_table)

        logger.info(f"✅ Saved: {filename}")
    except Exception as e:
        logger.error(f"❌ Error in statistical test for {filename}: {e}")

#******************************************************
# Funzione principale: mantiene il flusso originale dello script.
# 1) legge i conteggi raw
# 2) salva/merge i CSV
# 3) crea barplot e boxplot raw
# 4) crea normalizzazione e boxplot normalized
def main(data):
    print("\n############################# DESCRIPTIVE ANALYSIS #############################\n")
    logger.info("Start overview process...")
    path_merge_folder=os.path.join(data["Paths"]["output_folder"],"Merged_clean")
    logger.info(f"Merge data found in {path_merge_folder}")
    if not os.listdir(path_merge_folder):
        logger.critical(f"{path_merge_folder} is an empty directory")
    path_output_results=os.path.join(data["Paths"]["output_folder"],"Descriptive")
    logger.info(f"Descriptive output will be stored in {path_output_results}")
    dict_pheno_interested=data["Phenotypes"]["pheno_list"]
    logger.info(f"The following phenotypes will be considered: {dict_pheno_interested}")
    dict_raw_count=raw_count_cells(path_merge_folder,dict_pheno_interested)
    groups=list(dict_raw_count.keys())
    logger.info(f"{len(groups)} group(s) found: {groups}")
    create_output_dir(path_output_results,groups)
    # Statistics options are read from the config.
    # Existing keys: Stats -> p_adj and Stats -> sample_type.
    # New optional keys: p_adj_threshold, save_stats_csv, show_pvalues_on_boxplot.
    p_adjust, p_threshold, save_stats_csv, show_pvalues_on_boxplot, test = _get_stats_options(data)

    # CREATE CSV and merge
    phenos_list = list(dict_pheno_interested)
    create_and_merge_counts(path_output_results, dict_raw_count, phenos_list)


    # Plot e statistiche sui conteggi raw.
    if data["Descriptive"]["raw"]:
        bar_plot(path_output_results, dict_raw_count, "Raw")
        if len(groups)>=2:
            create_output_box_plot_dir(path_output_results)
            df_raw=prepare_data_box_plot(dict_raw_count)
            create_comparison_box_plot(path_output_results, df_raw, p_adjust, test, "Raw",
                                       p_threshold, save_stats_csv, show_pvalues_on_boxplot)
        else:
            logger.critical("Found only one group - impossible to continue with groups comparison and box plot figure")

    # Plot e statistiche sui conteggi normalizzati.
    if data["Descriptive"]["normalized"]:
        mean_group=calculate_mean_group_cells(dict_raw_count)
        norm_count=normalized_count_cells(dict_raw_count, mean_group)
        create_summury_file(path_output_results, norm_count, "Norm")
        merge_counts(path_output_results, "Norm_count_*.csv", "merged_norm.csv", phenos_list)
        bar_plot(path_output_results, norm_count, "Normalized")
        if len(groups)>=2:
            mean_group_all=calculate_mean_total_groups(dict_raw_count)
            norm_count_all_groups=normalized_count_on_all_groups(dict_raw_count, mean_group_all)
            create_norm_all_file(path_output_results, norm_count_all_groups)
            merge_counts(path_output_results, "all_norm_count_*.csv", "merged_all_norm.csv", phenos_list)
            create_output_box_plot_dir(path_output_results)
            df_norm=prepare_data_box_plot(norm_count_all_groups)
            create_comparison_box_plot(path_output_results, df_norm, p_adjust, test, "Normalized",
                                           p_threshold, save_stats_csv, show_pvalues_on_boxplot)
        else:
            logger.critical("Found only one group-impossible to continue with group comparison")

    logger.info("End descriptive analysis!\n")
    return()

if __name__=="__main__":
    main()
