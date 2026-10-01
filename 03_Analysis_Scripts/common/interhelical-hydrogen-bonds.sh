#!/bin/bash
# ============================================================
# Combined workflow for interhelical hydrogen bond analysis
# Steps:
#   1. Remove PBC jumps (nojump)
#   2. Create index groups for helices A-D
#   3. Compute interhelical H-bonds for six helix-helix pairs
#   4. Sum the six pairs and report total interhelical H-bonds
# ============================================================

# ========== User inputs ==========
TPR="md_500ns.tpr"
XTC_RAW="md.xtc"
XTC="md_nojump.xtc"

# Helix residue ranges (adjust to your system)
RES_A="2-27"
RES_B="32-56"
RES_C="59-76"
RES_D="80-99"

# Trajectory time range (ns) and output interval for hbond
B_TIME=0
E_TIME=500
DT=0.05

# Group numbers for helices (leave empty to auto-detect from index.ndx)
gA=""
gB=""
gC=""
gD=""
# =================================

LOGFILE="hbond_analysis.log"
> "$LOGFILE"

log() {
    echo "$1" | tee -a "$LOGFILE"
}

log "=========================================="
log "Interhelical hydrogen bond analysis"
log "Start: $(date)"
log "=========================================="

# ---------- Step 1: PBC nojump ----------
if [ ! -f "$XTC" ]; then
    log "Step 1: Removing PBC jumps -> $XTC"
    echo "0" | gmx trjconv -s "$TPR" -f "$XTC_RAW" -o "$XTC" -pbc nojump 2>&1 | tee -a "$LOGFILE"
else
    log "Step 1: $XTC already exists, skipping."
fi

# Check files
for f in "$TPR" "$XTC"; do
    if [ ! -f "$f" ]; then
        log "ERROR: File $f not found."
        exit 1
    fi
done

# ---------- Step 2: Create index file ----------
log "Step 2: Creating index.ndx with helix groups..."
rm -f index.ndx

echo -e "r $RES_A\nr $RES_B\nr $RES_C\nr $RES_D\nq" | \
gmx make_ndx -f "$TPR" -o index.ndx 2>&1 | tee -a "$LOGFILE"

if [ ! -f index.ndx ]; then
    log "ERROR: index.ndx creation failed."
    exit 1
fi

log "Groups in index.ndx:"
cat index.ndx | grep -E "^\[.*\]$" | nl -w2 -s': ' | tee -a "$LOGFILE"

# ---------- Auto-detect group numbers ----------
get_group_number() {
    local target="$1"
    awk -v target="$target" '
        /^\[/ {
            name=$0
            gsub(/^\[[[:space:]]*/, "", name)
            gsub(/[[:space:]]*\]$/, "", name)
            if (name == target) {
                print count
                exit
            }
            count++
        }
        BEGIN { count=0 }
    ' index.ndx
}

if [ -z "$gA" ]; then
    gA=$(get_group_number "r_${RES_A}")
fi
if [ -z "$gB" ]; then
    gB=$(get_group_number "r_${RES_B}")
fi
if [ -z "$gC" ]; then
    gC=$(get_group_number "r_${RES_C}")
fi
if [ -z "$gD" ]; then
    gD=$(get_group_number "r_${RES_D}")
fi

if [ -z "$gA" ] || [ -z "$gB" ] || [ -z "$gC" ] || [ -z "$gD" ]; then
    log "ERROR: Could not auto-detect all helix group numbers."
    log "Please check group names in index.ndx and set gA, gB, gC, gD manually."
    exit 1
fi

log "Using group numbers: A=$gA, B=$gB, C=$gC, D=$gD"

# ---------- Step 3: Compute six interhelical H-bond pairs ----------
PAIRS=("$gA $gB" "$gA $gC" "$gA $gD" "$gB $gC" "$gB $gD" "$gC $gD")
PAIR_NAMES=("AB" "AC" "AD" "BC" "BD" "CD")

for i in "${!PAIRS[@]}"; do
    pair=${PAIRS[$i]}
    name=${PAIR_NAMES[$i]}
    set -- $pair
    ref=$1
    target=$2
    outfile="pair_${name}.xvg"

    log "Computing H-bonds for helix pair $name (groups $ref vs $target)..."
    gmx hbond -f "$XTC" -s "$TPR" -n index.ndx \
        -r "$ref" -t "$target" \
        -num "$outfile" -tu ns -b "$B_TIME" -e "$E_TIME" -dt "$DT" 2>&1 | tee -a "$LOGFILE"

    if [ ! -f "$outfile" ]; then
        log "ERROR: $outfile was not generated."
        exit 1
    fi
    log "$outfile generated ✓"
done

# ---------- Step 4: Sum and statistics ----------
log "Summing six pairs and computing statistics..."
python3 << 'EOF' 2>&1 | tee -a "$LOGFILE"
import numpy as np

files = ["pair_AB.xvg", "pair_AC.xvg", "pair_AD.xvg",
         "pair_BC.xvg", "pair_BD.xvg", "pair_CD.xvg"]

time = None
total = None

for fname in files:
    data = np.loadtxt(fname, comments=['#', '@'])
    if time is None:
        time = data[:, 0]
        total = np.zeros(len(time))
    total += data[:, 1]

np.savetxt("total_cross_hbond.dat",
           np.column_stack((time, total)),
           fmt="%.6f %.0f",
           header="Time(ns)  Total_Interhelix_HBonds",
           comments='')

mean_val = np.mean(total)
std_val = np.std(total, ddof=1)
print("\n===== Final statistics =====")
print(f"Number of frames: {len(total)}")
print(f"Average interhelical H-bonds per frame: {mean_val:.2f} ± {std_val:.2f}")
print("Detailed data saved to total_cross_hbond.dat")
EOF

log "=========================================="
log "End: $(date)"
log "All done."
log "=========================================="