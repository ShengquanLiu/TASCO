import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.dccm import DynamicalCrossCorrelation

# =============================================
# USER CONFIGURATION
# =============================================

# True  -> compute DCCM from trajectories and save .npy files
# False -> skip computation, only plot existing .npy files
COMPUTE_DCCM = True

# True  -> plot the filtered matrix (only |i-j|>10 and |rho|>0.5 pairs kept)
# False -> plot the full DCCM matrix
PLOT_FILTERED_MATRIX = True

# Three replicates for each binder: list of (tpr, xtc) pairs.
# Replace with your actual file names.
TRAJECTORIES = {
    "IL7Ra-4H-01": [
        ("md_4H01_rep1.tpr", "md_4H01_rep1.xtc"),
        ("md_4H01_rep2.tpr", "md_4H01_rep2.xtc"),
        ("md_4H01_rep3.tpr", "md_4H01_rep3.xtc"),
    ],
    "IL7Ra-3H-01": [
        ("md_3H01_rep1.tpr", "md_3H01_rep1.xtc"),
        ("md_3H01_rep2.tpr", "md_3H01_rep2.xtc"),
        ("md_3H01_rep3.tpr", "md_3H01_rep3.xtc"),
    ],
    "IL7Rα-55": [
        ("md_a55_rep1.tpr", "md_a55_rep1.xtc"),
        ("md_a55_rep2.tpr", "md_a55_rep2.xtc"),
        ("md_a55_rep3.tpr", "md_a55_rep3.xtc"),
    ],
}

# Display names shown as row labels
DISPLAY_NAMES = {
    "IL7Ra-4H-01": "4H-01",
    "IL7Ra-3H-01": "3H-01",
    "IL7Rα-55":    "α-55",
}

# ---- Method-defined thresholds ----
SEQ_DIST_THRESHOLD = 10       # keep only pairs with |i - j| > 10
CORR_THRESHOLD = 0.5          # keep only pairs with |rho| > 0.5

# ---- Plotting parameters ----
VMIN, VMAX = -1, 1
CMAP = "RdBu_r"               # -1 -> blue, 0 -> white, +1 -> red
OUTPUT_TIF = "Figure_2F_all_replicates.tif"
OUTPUT_PDF = "Figure_2F_all_replicates.pdf"

# =============================================
# STEP 1: COMPUTE DCCM FOR EACH REPLICATE
# =============================================

def compute_dccm(tpr_file, xtc_file):
    """Load a trajectory, align it, and compute the DCCM using MDAnalysis."""
    u = mda.Universe(tpr_file, xtc_file)
    ca = u.select_atoms('protein and name CA')
    print(f"    Aligning trajectory (CA atoms: {len(ca)}) ...")
    align.AlignTraj(u, u, select='protein and name CA', in_memory=True).run()

    print(f"    Computing DCCM ...")
    dccm = DynamicalCrossCorrelation(u, select='name CA').run()
    return dccm.results.corr


def apply_method_filter(matrix,
                        seq_thresh=SEQ_DIST_THRESHOLD,
                        corr_thresh=CORR_THRESHOLD):
    """
    Apply Method-style filtering:
      - keep only pairs with |i - j| > seq_thresh
      - keep only pairs with |rho| > corr_thresh
    Everything else (including diagonal) is set to 0.

    Returns:
      filtered_matrix : same shape as input, non-qualifying entries zeroed
      strong_pairs    : list of dicts describing each retained pair
    """
    filtered = np.zeros_like(matrix)
    strong_pairs = []

    rows, cols = np.triu_indices_from(matrix, k=1)
    for i, j in zip(rows, cols):
        if abs(i - j) > seq_thresh and abs(matrix[i, j]) > corr_thresh:
            filtered[i, j] = matrix[i, j]
            filtered[j, i] = matrix[j, i]     # keep symmetry
            strong_pairs.append({
                "Residue_i": int(i),
                "Residue_j": int(j),
                "Correlation": float(matrix[i, j]),
                "Seq_Distance": int(abs(i - j))
            })

    np.fill_diagonal(filtered, 0)
    return filtered, strong_pairs


if COMPUTE_DCCM:
    for protein, reps in TRAJECTORIES.items():
        print(f"\n{'=' * 60}")
        print(f"Processing: {protein}")
        print(f"{'=' * 60}")

        for i, (tpr, xtc) in enumerate(reps, start=1):
            print(f"  Replicate {i}: {xtc}")
            mat = compute_dccm(tpr, xtc)

            # ---- Full DCCM ----
            full = mat.copy()
            np.fill_diagonal(full, 0)
            np.save(f"{protein}_rep{i}_full.npy", full)
            print(f"    ✅ Saved full DCCM: {protein}_rep{i}_full.npy")

            # ---- Filtered DCCM ----
            filtered, strong_pairs = apply_method_filter(full)
            np.save(f"{protein}_rep{i}_filtered.npy", filtered)
            print(f"    ✅ Saved filtered DCCM: {protein}_rep{i}_filtered.npy")

            # ---- Strong-pair list ----
            pd.DataFrame(
                strong_pairs,
                columns=["Residue_i", "Residue_j", "Correlation", "Seq_Distance"]
            ).to_csv(f"{protein}_rep{i}_Strong_Pairs.csv", index=False)
            print(f"    ✅ Strong pairs (|rho|>{CORR_THRESHOLD}, "
                  f"dist>{SEQ_DIST_THRESHOLD}): {len(strong_pairs)}")
            print(f"       Saved: {protein}_rep{i}_Strong_Pairs.csv")

# =============================================
# STEP 2: PLOT 3 x 3 GRID (rows = protein, cols = replicate)
# =============================================

plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

proteins = list(DISPLAY_NAMES.keys())
n_rows = len(proteins)
n_cols = 3

fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 18))

for r, protein in enumerate(proteins):
    display = DISPLAY_NAMES[protein]

    for c in range(n_cols):
        ax = axes[r, c]
        rep = c + 1
        suffix = "filtered" if PLOT_FILTERED_MATRIX else "full"
        filepath = f"{protein}_rep{rep}_{suffix}.npy"

        if not os.path.exists(filepath):
            ax.text(0.5, 0.5, f"File not found:\n{filepath}",
                    ha='center', va='center')
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        matrix = np.load(filepath)
        n_res = matrix.shape[0]

        sns.heatmap(matrix, ax=ax, vmin=VMIN, vmax=VMAX, cmap=CMAP,
                    center=0, square=True,
                    xticklabels=20 if n_res > 80 else 10,
                    yticklabels=20 if n_res > 80 else 10,
                    cbar=False)

        # Column title (replicate label) on top row only
        if r == 0:
            ax.set_title(f"rep{c + 1}", fontsize=18, pad=10)

        # Row label (protein name) on first column only
        if c == 0:
            ax.set_ylabel(f"{display}\nResidue No.", fontsize=16)
        else:
            ax.set_ylabel("Residue No.", fontsize=16)

        ax.set_xlabel("Residue No.", fontsize=16)

        # Tick positions and labels (same as Figure2D-dccm-plot.py)
        step = 20 if n_res > 80 else 10
        tick_positions = np.arange(0, n_res, step) + 0.5
        tick_labels = np.arange(0, n_res, step)
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=15)
        ax.set_yticklabels(tick_labels, fontsize=15)

        # NOTE: No cleavage-site markers are drawn here.
        # The heatmap shows only the DCCM matrix values.

# Unified colorbar on the right
cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
norm = plt.Normalize(vmin=VMIN, vmax=VMAX)
sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label("Correlation (ρ)", fontsize=18)

plt.subplots_adjust(left=0.06, right=0.90, wspace=0.30, hspace=0.30, bottom=0.06)
plt.savefig(OUTPUT_TIF, dpi=600, format='tif', bbox_inches='tight')
plt.savefig(OUTPUT_PDF, format='pdf', bbox_inches='tight')
print(f"\n✅ Saved: {OUTPUT_TIF}")
print(f"✅ Saved: {OUTPUT_PDF}")
plt.show()