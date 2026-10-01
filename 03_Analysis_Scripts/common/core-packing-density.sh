# 1. Handle PBC (skip if already done)
gmx trjconv -f md.xtc -s md_500ns.tpr -o md_noPBC.xtc -pbc mol -center
# Select Protein, then System

# 2. Build hydrophobic core index (using the residue numbers you determined)
gmx sasa -f md_noPBC.xtc -s md_500ns.tpr -surface Protein -output Sidechain -or residue_area.xvg
# Residues with average side-chain SASA < 10 Å² are defined as core residues (a widely cited classical threshold).
awk 'NR>1 && $1 !~ /^#/ && $1 !~ /^@/ && $4 < 0.10 {print $1}' residue_area.xvg > core_residues.txt

gmx make_ndx -f md_500ns.tpr -o hydrophobic_core.ndx

r 7 18 25 30 39 42 43 46 50 53 61 62 65 68 72 75 83 87 90 91 94 97
name <group name> Hydrophobic_Core
q

# 3. Extract side-chain atom coordinates and calculate V_vdW
gmx trjconv -s md_500ns.tpr -f md_noPBC.xtc -n hydrophobic_core.ndx -o core_sidechain.gro -dump 0
# Select Hydrophobic_Core

V_vdW=$(awk 'NR>2 && NF>=6 {ele=substr($2,1,1); if(ele=="C"||ele=="N"||ele=="O"||ele=="S"||ele=="H") c[ele]++} END {pi=3.1415926535; v=c["C"]*(4/3*pi*0.170^3)+c["N"]*(4/3*pi*0.155^3)+c["O"]*(4/3*pi*0.140^3)+c["S"]*(4/3*pi*0.180^3)+c["H"]*(4/3*pi*0.120^3); printf "%.6f", v}' core_sidechain.gro)
echo "V_vdW = $V_vdW nm^3"

# 4. Calculate the total volume V_total over time
gmx sasa -f md_noPBC.xtc -s md_500ns.tpr -n hydrophobic_core.ndx -surface Hydrophobic_Core -tv volume.xvg -o area.xvg
# Select Hydrophobic_Core

# 5. Divide frame by frame to obtain the packing fraction trajectory
awk -v vdw=$V_vdW 'NR>1 && $1 !~ /^#/ && $1 !~ /^@/ {printf "%.4f %.6f\n", $1, vdw/$2}' volume.xvg > packing_density.xvg
echo "Done: packing_density.xvg"