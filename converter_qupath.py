import pandas as pd
import os

# ==============================
# CONFIGURAZIONE
# ==============================

input_folder = "/Users/elettramercante/Desktop/Integrazione_Qupath_ALOA/cut"
output_folder = "/Users/elettramercante/Desktop/ALOA/data_test/raw_data/TEST_QUPATH"
os.makedirs(output_folder, exist_ok=True)

markers = ["PANCK", "CD4", "CD8", "CD20", "CD66b", "CD68"]

output_path = os.path.join(
    output_folder,
    "TEST_QUPATH_[00000,00000]_cell_seg_data.txt"
)

# ==============================
# FUNZIONE LETTURA QUPATH
# ==============================

def load_qupath(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except:
        df = pd.read_csv(path, sep=",")

    df = df.rename(columns={
        "Object ID": "Cell ID",
        "Centroid X µm": "Cell X Position",
        "Centroid Y µm": "Cell Y Position",
        "ROI": "Process Region ID"
    })

    if "Classification" not in df.columns:
        raise RuntimeError(f"❌ Colonna 'Classification' mancante in {path}")

    return df


# ==============================
# FILE BASE (struttura celle)
# ==============================

base_marker = markers[0]
base_file = os.path.join(input_folder, f"{base_marker.lower()}_detection_cut.txt")

df_base = load_qupath(base_file)

df_final = pd.DataFrame()
df_final["Cell ID"] = df_base["Cell ID"]
df_final["Cell X Position"] = df_base["Cell X Position"]
df_final["Cell Y Position"] = df_base["Cell Y Position"]
df_final["Process Region ID"] = df_base.get("Process Region ID", "#N/A")

df_final["Path"] = input_folder
# ⚠️ NOMI ESATTI come ALOA li vuole
df_final["Sample.Name"] = "TEST_QUPATH_[00000,00000]"
df_final["Slide.ID"] = "TEST_QUPATH_[00000,00000]"

df_final["Total cells"] = 1
df_final["Cell Density (per square mm)"] = "#N/A"

# ==============================
# PHENOTYPES MULTI-MARKER
# ==============================

for marker in markers:
    file_path = os.path.join(input_folder, f"{marker.lower()}_detection_cut.txt")
    df_m = load_qupath(file_path)

    df_final[f"Phenotype-{marker}"] = df_m["Classification"].astype(str).apply(
        lambda x: f"{marker}+" if x.strip().upper() == marker else "OTHER"
    )

# ==============================
# ORDINE COLONNE (ALOA)
# ==============================

ordered_cols = (
    ["Path", "Sample.Name", "Slide.ID"]
    + [f"Phenotype-{m}" for m in markers]
    + [
        "Cell ID",
        "Total cells",
        "Cell Density (per square mm)",
        "Cell X Position",
        "Cell Y Position",
        "Process Region ID",
    ]
)

df_final = df_final[ordered_cols]

# ==============================
# OUTPUT
# ==============================

df_final.to_csv(output_path, sep="\t", index=False)

print("✔️ Conversione QuPath → ALOA completata")
print(f"📄 File salvato in:\n{output_path}")

