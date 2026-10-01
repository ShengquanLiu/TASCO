import mdtraj as md
import numpy as np
from tqdm import tqdm

# ================= Configuration =================
traj_file = 'md_center.xtc'
top_file  = 'protein_only.gro'
out_prefix = 'rep1'
probe_radius = 0.14   # nm, i.e. 1.4 Å
# =================================================

# ---------- Miller et al. 1987 maximum side-chain SASA (Å²) ----------
max_sidechain_sasa_miller = {
    'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0,
    'CYS': 140.0, 'GLN': 189.0, 'GLU': 183.0, 'GLY':  85.0,
    'HIS': 194.0, 'ILE': 182.0, 'LEU': 180.0, 'LYS': 211.0,
    'MET': 204.0, 'PHE': 218.0, 'PRO': 143.0, 'SER': 122.0,
    'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0,
    # Non-standard residue names mapped to the corresponding Miller values
    'HID': 194.0, 'HIE': 194.0, 'HIP': 194.0,
    'CYX': 140.0,
    'ASH': 151.0, 'GLH': 183.0, 'LYN': 211.0,
    'ASPH': 151.0, 'GLUH': 183.0,
}

# ---------- 1. Load trajectory and remove hydrogens ----------
traj = md.load(traj_file, top=top_file)
traj = traj.remove_hydrogens()   # consistent with the heavy-atom reference of Miller 1987
n_frames   = traj.n_frames
n_residues = traj.n_residues
top = traj.topology
print(f"Trajectory: {n_frames} frames, {n_residues} residues, {traj.n_atoms} atoms (hydrogens removed)")

# ---------- 2. Select side-chain heavy atoms ----------
# Backbone heavy-atom names in the GROMACS amber force field: N, CA, C, O, OXT
backbone_names = {'N', 'CA', 'C', 'O', 'OXT'}

res_names = []
res_sidechain_atoms = []

for ri in range(n_residues):
    res = top.residue(ri)
    res_names.append(res.name)
    atoms = np.array([a.index for a in res.atoms if a.name not in backbone_names])
    res_sidechain_atoms.append(atoms)

res_names = np.array(res_names)

# Check for residues without a Miller reference value
missing = set(res_names) - set(max_sidechain_sasa_miller.keys())
if missing:
    print(f"Warning: the following residue types have no Miller reference value and will be set to NaN: {missing}")

# ---------- 3. Compute side-chain SASA per frame ----------
sidechain_sasa_A2 = np.zeros((n_frames, n_residues), dtype=np.float32)

print("Computing side-chain SASA...")
for i in tqdm(range(n_frames)):
    sasa_atom = md.shrake_rupley(traj[i], mode='atom', probe_radius=probe_radius).flatten()
    for ri in range(n_residues):
        atoms = res_sidechain_atoms[ri]
        if len(atoms) > 0:
            sidechain_sasa_A2[i, ri] = sasa_atom[atoms].sum() * 100.0  # nm² -> Å²
        else:
            sidechain_sasa_A2[i, ri] = 0.0

# ---------- 4. Save raw side-chain SASA, residue names and indices ----------
np.save(f'{out_prefix}_sidechain_sasa_A2_all_residues.npy', sidechain_sasa_A2)
np.save(f'{out_prefix}_residue_names_all.npy', res_names)
np.save(f'{out_prefix}_residue_indices_all.npy', np.arange(n_residues))

print(f"Saved: {out_prefix}_sidechain_sasa_A2_all_residues.npy")
print(f"   Shape: {sidechain_sasa_A2.shape}  (n_frames, n_residues)")

# ---------- 5. Normalize with the Miller dictionary to obtain rSASA ----------
rsasa_miller = np.full_like(sidechain_sasa_A2, np.nan, dtype=np.float32)

for ri in range(n_residues):
    name = str(res_names[ri])
    if name in max_sidechain_sasa_miller and max_sidechain_sasa_miller[name] > 0:
        rsasa_miller[:, ri] = sidechain_sasa_A2[:, ri] / max_sidechain_sasa_miller[name]
    else:
        rsasa_miller[:, ri] = np.nan

# ---------- 6. Save Miller-normalized rSASA ----------
np.save(f'{out_prefix}_rsasa_Miller_all_residues.npy', rsasa_miller)
print(f"Saved: {out_prefix}_rsasa_Miller_all_residues.npy")
print(f"   Shape: {rsasa_miller.shape}  (n_frames, n_residues)")

# ---------- 7. Optional: quick exposure statistics (percentage of frames with rSASA > 0.3) ----------
# Note: this assumes equally spaced frames; if the frame spacing is uneven, weight by time instead.
exposed_fraction = np.nanmean(rsasa_miller > 0.3, axis=0) * 100
# Target sites can be selected later based on residue indices, or saved to file.
np.save(f'{out_prefix}_exposed_fraction_all_residues.npy', exposed_fraction)
print(f"Saved: {out_prefix}_exposed_fraction_all_residues.npy")
print("   Exposure percentages (%) were computed with Miller normalization and can be filtered for target residues later.")