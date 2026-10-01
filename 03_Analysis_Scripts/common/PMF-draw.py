import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re

# Automatically get the folder path where the current script is located
current_dir = os.path.dirname(os.path.abspath(__file__))

# Example: For IL7Ra-4H-01 site14, the PMF output files follow this naming pattern:
#          pmf-site14_final_t*.xvg (e.g., pmf-site14_final_t50000.xvg, pmf-site14_final_t49800.xvg, ...)

file_pattern = os.path.join(current_dir, "pmf-site14_final_t*.xvg")
save_path = os.path.join(current_dir, "PMF_Site14_conv.png")

def read_xvg(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line[0] in ('#', '@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    pass
    return np.array(data)

files = sorted(glob.glob(file_pattern))
print(f"Number of files found: {len(files)}")

if not files:
    print("No files found. Please check the filename or path!")
else:
    # Keep canvas size unchanged
    plt.figure(figsize=(14, 10), dpi=600)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

    for i, f in enumerate(files):
        data = read_xvg(f)
        time_match = re.search(r't(\d+)', os.path.basename(f))
        if time_match:
            time_val = float(time_match.group(1)) / 1000.0
            label = f"{time_val:.1f} ns"
        else:
            label = f"Run {i+1}"

        plt.plot(data[:, 0], data[:, 1], linewidth=4, label=label, color=colors[i % len(colors)], alpha=0.9)

    # ---- Axis labels ----
    # Keep font size 42 bold; increase labelpad to prevent covering the ticks
    plt.xlabel("Reaction Coordinate (nm)", fontsize=42, fontweight='bold', labelpad=40)
    plt.ylabel("PMF (kJ/mol)", fontsize=42, fontweight='bold', labelpad=40)

    # ---- Tick values: increased from 30 to 36 ----
    plt.tick_params(axis='both', which='major', labelsize=36, width=4, length=10, direction='out')

    # Force tick numbers to be bold
    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Make border thicker
    for spine in ax.spines.values():
        spine.set_linewidth(4)

    # Keep legend at size 30 bold to avoid overpowering
    plt.legend(fontsize=30, loc='lower left', frameon=True, framealpha=0.9, edgecolor='black', prop={'weight': 'bold', 'size': 30})

    # Force uniform Y-axis
    plt.ylim(0, 50)

    # Make grid lines thicker
    plt.grid(True, linestyle='--', alpha=0.6, linewidth=2)

    plt.tight_layout()
    plt.savefig(save_path, dpi=600, format='png', bbox_inches='tight')
    print(f"Image saved to: {save_path}")

    try:
        plt.show()
    except:
        print("Unable to open a window automatically. Please check the saved PNG image in the folder.")