import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import pandas as pd

# =============================================
# USER CONFIGURATION
# =============================================

PROTEINS = {
    "IL7Ra-4H-01": {
        "input_file":     "IL7Ra-4H-01.csv",
        "output_file":    "Figure_S7-IL7Ra-4H-01.tiff",
        "total_residues": 100,
        "leu_residues":   [25, 42, 49, 58, 64, 65, 71, 84, 86, 87, 91, 93, 94],
        "met_residues":   [14],
        "tyr_residues":   [72],
        "phe_residues":   [],          # no Phe arrows for this protein
    },
    "IL7Rα-55": {
        "input_file":     "IL7Rα-55.csv",
        "output_file":    "Figure_S8-IL7Rα-55.tiff",
        "total_residues": 65,
        "leu_residues":   [3, 53],
        "met_residues":   [52],
        "tyr_residues":   [],
        "phe_residues":   [62],
    },
    "IL7Ra-3H-01": {
        "input_file":     "IL7Ra-3H-01.csv",
        "output_file":    "Figure_S8-IL7Ra-3H-01.tiff",
        "total_residues": 71,
        "leu_residues":   [8, 23, 32, 39, 43, 50, 65, 69],
        "met_residues":   [42],
        "tyr_residues":   [22, 31, 46, 49],
        "phe_residues":   [],
    },
}

# Arrow colors
COLOR_L = '#000000'   # Leu -> black
COLOR_M = '#8A2BE2'   # Met -> purple
COLOR_Y = '#2ca02c'   # Tyr -> green
COLOR_F = '#FF8C00'   # Phe -> dark orange

KR_FONTSIZE     = 11
LEGEND_FONTSIZE = 13

# =============================================
# HELPERS
# =============================================

dssp_map = {
    'H': 1, 'E': 2, 'T': 3,
    'G': 4, 'I': 4, 'B': 4, 'S': 4, 'C': 4,
    '~': 4, ' ': 4, '':  4,
}

def to_num(x):
    if pd.isna(x):
        return 4
    s = str(x).strip().upper()
    if s in dssp_map:
        return dssp_map[s]
    try:
        v = float(s)
        if 1 <= v <= 4:
            return int(v)
    except ValueError:
        pass
    return 4


def res_to_y(n):
    return n - 0.5


def add_arrows(ax, residues, label, color, total_residues):
    valid = [r for r in residues if 1 <= r <= total_residues]
    for i, res in enumerate(valid):
        y = res_to_y(res)
        tier = i % 3
        x_text = 505 + tier * 22
        ax.annotate(
            f'{label}{res}',
            xy=(500, y),
            xytext=(x_text, y),
            ha='left', va='center',
            fontsize=KR_FONTSIZE,
            fontweight='bold', color=color,
            annotation_clip=False,
            arrowprops=dict(arrowstyle='-|>',
                            color=color, lw=1.3,
                            shrinkA=1, shrinkB=0)
        )


def plot_one_protein(cfg):
    input_file     = cfg["input_file"]
    output_file    = cfg["output_file"]
    total_residues = cfg["total_residues"]

    # ---- Read data ----
    try:
        df = pd.read_csv(input_file, header=None, dtype=str)
        if df.shape[1] < 2:
            raise ValueError
    except Exception:
        df = pd.read_csv(input_file, header=None, dtype=str,
                         sep=r'\s+', engine='python')

    print(f"\n{'=' * 60}")
    print(f"Protein: {input_file}")
    print(f"CSV shape: {df.shape}")

    time_series = df.iloc[:, 0].astype(str).str.extract(r'([\d.]+)')[0]
    time = pd.to_numeric(time_series, errors='coerce').fillna(0).values.astype(float)

    data   = df.iloc[:, 1:].apply(lambda col: col.map(to_num)).values.astype(float)
    data_t = data.T

    print(f"Number of residue columns in CSV = {data_t.shape[0]}  (expected = {total_residues})")
    if data_t.shape[0] != total_residues:
        print("⚠️ Number of data columns does not match total_residues, please check!")

    vals, counts = np.unique(data_t, return_counts=True)
    print("State distribution:", dict(zip(vals.tolist(), counts.tolist())))

    # ---- Color mapping ----
    colors = ['#E31A1C', '#1F5FA8', '#FFD700', '#BFBFBF']
    cmap   = ListedColormap(colors)
    bounds = [0.5, 1.5, 2.5, 3.5, 4.5]
    norm   = BoundaryNorm(bounds, cmap.N)

    # ---- Plotting ----
    fig, ax = plt.subplots(figsize=(12, 6))

    im = ax.imshow(data_t, aspect='auto', cmap=cmap, norm=norm,
                   extent=[time.min(), time.max(), 0, total_residues],
                   origin='lower', interpolation='none')

    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_ylabel('Residue number', fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.set_xlim(0, 500)
    ax.set_ylim(0, total_residues)

    # ---- Leu / Met / Tyr / Phe arrow annotations ----
    add_arrows(ax, cfg["met_residues"], 'M', COLOR_M, total_residues)
    add_arrows(ax, cfg["leu_residues"], 'L', COLOR_L, total_residues)
    add_arrows(ax, cfg["tyr_residues"], 'Y', COLOR_Y, total_residues)
    add_arrows(ax, cfg["phe_residues"], 'F', COLOR_F, total_residues)

    # ---- Legend ----
    legend_elements = [
        Patch(facecolor='#E31A1C', label='α-Helix'),
        Patch(facecolor='#1F5FA8', label='β-Sheet'),
        Patch(facecolor='#FFD700', label='Turn'),
        Patch(facecolor='#BFBFBF', label='Coil'),
    ]
    ax.legend(handles=legend_elements,
              loc='upper left',
              bbox_to_anchor=(1.18, 1),
              prop={'size': LEGEND_FONTSIZE},
              frameon=False,
              handlelength=1.8,
              handleheight=1.1,
              borderaxespad=0.0,
              labelspacing=0.8)

    # ---- Save ----
    plt.tight_layout()
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.show()
    print(f"Figure saved as: {output_file}")


# =============================================
# RUN ALL
# =============================================

if __name__ == "__main__":
    for protein_name, cfg in PROTEINS.items():
        plot_one_protein(cfg)