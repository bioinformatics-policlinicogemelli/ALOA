#SPIEGAZIONE: script che genera 3 tipologie di heatmaps, in base alla scelta nel config

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

# righe sottili per celle non significative
mpl.rcParams["hatch.linewidth"] = 0.2


# -------------------
# Lettura CSV
# -------------------
def load_csv(pairwise_file: str) -> pd.DataFrame:
    """Carica il CSV in un DataFrame pandas."""
    return pd.read_csv(pairwise_file)


# -------------------
# Costruzione matrici heatmap
# -------------------
def build_heatmap(df: pd.DataFrame,
                  g1: str,
                  g2: str,
                  phenos_from,
                  phenos_to,
                  alpha: float = 0.05):

    hm1 = pd.DataFrame(np.nan, index=phenos_from, columns=phenos_to)
    hm2 = pd.DataFrame(np.nan, index=phenos_from, columns=phenos_to)

    sub = df[(df["GROUP_1"] == g1) & (df["GROUP_2"] == g2)]

    p_col = "P_VALUE_ADJ" if "P_VALUE_ADJ" in df.columns else "P_VALUE"

    for _, row in sub.iterrows():
        f = row["PHENO_FROM"]
        t = row["PHENO_TO"]
        pval = row.get(p_col, np.nan)

        if pd.notna(pval) and pval < alpha:
            hm1.loc[f, t] = row["MEDIAN_1"]
            hm2.loc[f, t] = row["MEDIAN_2"]
        else:
            hm1.loc[f, t] = np.nan
            hm2.loc[f, t] = np.nan

    return hm1, hm2


# -------------------
# Annotazioni numeriche con colore adattivo
# -------------------
def annotate_matrix(ax, hm: pd.DataFrame, im, fontsize: float):
    """
    Scrive il valore numerico dentro le celle non-NaN. Il colore del testo è adattato al colore della cella: testo nero su celle chiare, bianco su celle scure.
    """
    n_rows, n_cols = hm.shape

    for i in range(n_rows):
        for j in range(n_cols):
            val = hm.iloc[i, j]
            if np.isnan(val):
                continue

            rgba = im.cmap(im.norm(val))
            r, g, b, _ = rgba

            luminance = 0.299 * r + 0.587 * g + 0.114 * b
            text_color = "black" if luminance > 0.6 else "white"

            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                fontsize=fontsize,
                color=text_color,
            )


# -------------------
# Funzione per creare una figura heatmap
# -------------------
def make_heatmap_figure(hm1: pd.DataFrame,
                        hm2: pd.DataFrame,
                        phenos_from,
                        phenos_to,
                        g1: str,
                        g2: str,
                        mode: str,
                        out_path: str,
                        cell_size: float = 0.5,
                        do_annot: bool = True):
    """
    crea figura con due heatmap affiancate (sx: g1, dx: g2)
    mode:
    -adaptive / adaptive_annot --> scala adattiva (vmin/vamx dal dato)
    -fixed / fixed_annot --> scala fissa [-1, 0, 1]
    """

    # dimensione automatica in base a n_fenotipi
    n_rows, n_cols = hm1.shape

    # dimensione desiderata di UNA cella (in pollici)
    heatmap_width = cell_size * n_cols
    heatmap_height = cell_size * n_rows

    # margini in pollici
    left_margin = 1.8
    right_margin = 1.2
    top_margin = 2.0      # spazio per colorbar + titoli
    bottom_margin = 2.0   # spazio per etichette X
    gap_between = 1.0     # spazio fra le due heatmap

    fig_width = max(8, 2 * heatmap_width + gap_between + left_margin + right_margin)
    fig_height = max(6, heatmap_height + top_margin + bottom_margin)

    fig = plt.figure(figsize=(fig_width, fig_height))

    # margini in coordinate [0,1] della figura
    left_frac = left_margin / fig_width
    right_frac = 1 - right_margin / fig_width
    bottom_frac = bottom_margin / fig_height
    top_frac = 1 - top_margin / fig_height

    # gridspec: una riga, due colonne (le due heatmap)
    gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        wspace=0.1,
        left=left_frac,
        right=right_frac,
        bottom=bottom_frac,
        top=top_frac,
    )

    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1], sharey=ax_left)

    for ax in (ax_left, ax_right):
        ax.set_aspect("auto")

    # Colormap: negativo = rosso, positivo = blu, NaN = grigio
    cmap = plt.cm.RdBu.copy()
    cmap.set_bad(color="lightgrey")

    # estrai min e max (ignorando i NaN) per le versioni adattive
    data_vals = np.concatenate([hm1.values.ravel(), hm2.values.ravel()])
    data_vals = data_vals[~np.isnan(data_vals)]
    if data_vals.size == 0:
        print("Nessun valore valido per la heatmap, skip.")
        plt.close(fig)
        return

    vmin = float(data_vals.min())
    vmax = float(data_vals.max())
    if vmin == vmax:
        vmin -= 1.0
        vmax += 1.0

    # Scelta della normalizzazione
    if mode.startswith("adaptive"):
        # colormap adattata al range del dato
        norm = None
    else:
        # scala fissa [-1, 0, 1]
        norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)

    # Disegno delle due heatmap
    if norm is None:
        im1 = ax_left.imshow(hm1, cmap=cmap, vmin=vmin, vmax=vmax)
        im2 = ax_right.imshow(hm2, cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        im1 = ax_left.imshow(hm1, cmap=cmap, norm=norm)
        im2 = ax_right.imshow(hm2, cmap=cmap, norm=norm)

    # Titoli e assi
    ax_left.set_title(str(g1), fontsize=15, fontweight="bold")
    ax_right.set_title(str(g2), fontsize=15, fontweight="bold")

    ax_left.set_xticks(range(len(phenos_to)))
    ax_left.set_xticklabels(phenos_to, rotation=90)
    ax_left.set_yticks(range(len(phenos_from)))
    ax_left.set_yticklabels(phenos_from)

    ax_right.set_xticks(range(len(phenos_to)))
    ax_right.set_xticklabels(phenos_to, rotation=90)
    ax_right.tick_params(axis="y", which="both", left=False, labelleft=False)

    # Celle non significative: NaN → righette bianche
    nan_mask1 = np.isnan(hm1.values)
    nan_mask2 = np.isnan(hm2.values)

    for i in range(nan_mask1.shape[0]):
        for j in range(nan_mask1.shape[1]):
            if nan_mask1[i, j]:
                rect = Rectangle(
                    (j - 0.5, i - 0.5),
                    1, 1,
                    facecolor="none",
                    hatch="///",
                    edgecolor="white",
                    linewidth=0,
                    zorder=3,
                )
                ax_left.add_patch(rect)

    for i in range(nan_mask2.shape[0]):
        for j in range(nan_mask2.shape[1]):
            if nan_mask2[i, j]:
                rect = Rectangle(
                    (j - 0.5, i - 0.5),
                    1, 1,
                    facecolor="none",
                    hatch="///",
                    edgecolor="white",
                    linewidth=0,
                    zorder=3,
                )
                ax_right.add_patch(rect)

    # Font size adattivo (in base alla dimensione della cella)
    fontsize = cell_size * 72 * 0.45  # ~45% dell'altezza cella
    fontsize = max(5, min(11, fontsize))

    # annotazioni solo se richiesto
    if do_annot:
        annotate_matrix(ax_left, hm1, im1, fontsize)
        annotate_matrix(ax_right, hm2, im2, fontsize)

    # Colorbar sopra, centrata e più corta
    fig.canvas.draw()  # aggiornamento posizioni assi

    bbox_left = ax_left.get_position()
    bbox_right = ax_right.get_position()

    x0 = bbox_left.x0
    x1 = bbox_right.x1
    full_width = x1 - x0

    # fattore di accorciamento: 0.5 = 50% della larghezza fra le due heatmap
    factor = 0.5
    cbar_width = full_width * factor
    center = (x0 + x1) / 2
    cbar_x0 = center - cbar_width / 2

    # posizione verticale leggermente sopra le heatmap
    cbar_height = 0.02
    gap = 0.15          # spazio tra heatmap e colorbar
    cbar_y0 = bbox_left.y1 + gap

    cbar_ax = fig.add_axes([cbar_x0, cbar_y0, cbar_width, cbar_height])
    cbar = fig.colorbar(im1, cax=cbar_ax, orientation="horizontal")

    # Per la versione a scala fissa, imposta tick fissi
    if mode.startswith("fixed"):
        cbar.set_ticks([-1.0, 0.0, 1.0])
        cbar.set_ticklabels(["-1", "0", "1"])

    # Salvataggio figura
    fig.savefig(out_path, dpi=600, bbox_inches="tight", pad_inches=0.30)
    print(f"Salvata: {out_path}")
    plt.close(fig)


# -------------------
# FUNZIONE PRINCIPALE 
# -------------------
def plot_all_heatmaps(data: dict,  
                      alpha: float = 0.05,
                      cell_size: float = 0.5
                      ):
    
    """
    Usa ALL_PAIRWISE_RESULTS.csv (cercato in Distance_Statistical/pairwise_pvalues) e il config 'data' per
    generare le heatmap in:

      <.../Distance_Statistical>/heatmaps_nuove/<GROUP1>_vs_<GROUP2>/

    Le scelte sono prese da:
      - data["Distance"]["heatmap_annot"]       (True/False)
      - data["Distance"]["heatmap_fixed_scale"] (True/False)
    """
    
    # leggo scelte dal config
    dist_cfg = data.get("Distance", {})
    heatmap_annot = bool(dist_cfg.get("heatmap_annot", True))
    heatmap_fixed_scale = bool(dist_cfg.get("heatmap_fixed_scale", True))

    # path del file ALL_PAIRWISE_RESULTS.csv
    output_folder = data["Paths"]["output_folder"]
    distance_stat_dir = os.path.join(output_folder, "Distance_Statistical")
    pairwise_dir = os.path.join(distance_stat_dir, "pairwise_pvalues")
    pairwise_file = os.path.join(pairwise_dir, "ALL_PAIRWISE_RESULTS.csv")
  
    if not os.path.isfile(pairwise_file):
        print(f"File pairwise non trovato: {pairwise_file}. Nessuna heatmap generata.")
        return

    df = load_csv(pairwise_file)

    for col in ["P_VALUE", "P_VALUE_ADJ", "MEDIAN_1", "MEDIAN_2"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # ordine fenotipi
    phenos = sorted(set(df["PHENO_FROM"]).union(df["PHENO_TO"]))
    if not phenos:
        print("Nessun fenotipo trovato nel CSV, controllo colonne PHENO_FROM / PHENO_TO.")
        return

    # coppie di gruppi presenti nel ALL_PAIRWISE_RESULTS.csv
    pairs_df = df[["GROUP_1", "GROUP_2"]].drop_duplicates().dropna()
    pairs = list(pairs_df.itertuples(index=False, name=None))
    if not pairs:
        print("Nessuna coppia GROUP_1 / GROUP_2 trovata nel CSV.")
        return

    base_name = os.path.splitext(os.path.basename(pairwise_file))[0]

    # cartella principale heatmaps SEMPRE dentro Distance_Statistical
    root_out_dir = os.path.join(distance_stat_dir, "heatmaps_nuove")
    os.makedirs(root_out_dir, exist_ok=True)

    print("Coppie di gruppi trovate:", pairs)
    print(f"Output in: {root_out_dir}\n")

    for g1, g2 in pairs:
        hm1, hm2 = build_heatmap(df, g1, g2, phenos, phenos, alpha=alpha)

        if np.isnan(hm1.values).all() and np.isnan(hm2.values).all():
            print(f"Nessun valore significativo per {g1} vs {g2}, skip.")
            continue

        pair_label = f"{g1}_vs_{g2}"
        pair_dir = os.path.join(root_out_dir, pair_label)
        os.makedirs(pair_dir, exist_ok=True)

        # 1) SEMPRE: colormap adattativa, senza numeri
        out_adapt = os.path.join(pair_dir, f"{base_name}_adaptive.png")
        make_heatmap_figure(
            hm1, hm2, phenos, phenos, g1, g2,
            mode="adaptive",
            out_path=out_adapt,
            cell_size=cell_size,
            do_annot=False,
        )

        # 2) SE heatmap_annot: adattiva + numeri
        if heatmap_annot:
            out_adapt_ann = os.path.join(pair_dir, f"{base_name}_adaptive_annot.png")
            make_heatmap_figure(
                hm1, hm2, phenos, phenos, g1, g2,
                mode="adaptive_annot",
                out_path=out_adapt_ann,
                cell_size=cell_size,
                do_annot=True,
            )

        # 3) SE heatmap_fixed_scale: colorbar fissa [-1,0,1]
        if heatmap_fixed_scale:
            if heatmap_annot:
                fname_fixed = f"{base_name}_fixed_annot.png"
                do_annot_fixed = True
            else:
                fname_fixed = f"{base_name}_fixed.png"
                do_annot_fixed = False

            out_fixed = os.path.join(pair_dir, fname_fixed)
            make_heatmap_figure(
                hm1, hm2, phenos, phenos, g1, g2,
                mode="fixed_annot",                # o "fixed": serve solo per usare la scala fissa
                out_path=out_fixed,
                cell_size=cell_size,
                do_annot=do_annot_fixed,
            )