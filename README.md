# Molecular Dynamics Data for IL7Rα Variants at pH 2 and pH 7

## Overview

This repository contains input files, key GROMACS outputs, analysis scripts, and processed data for the molecular dynamics (MD) simulations of IL7Rα variants under different pH conditions. Each system was simulated in three independent 500-ns replicates.

The variants included are:
- IL7Ra-4H-01 (pH 2 and pH 7)
- IL7Ra-3H-01 (pH 2 only)
- IL7Rα-55 (pH 2 only)
- IL7Ra-4H-RD40 (pH 2 and pH 7)
- IL7Ra-4H-RD30 (pH 2 and pH 7)

All simulations were performed with GROMACS 2024.6 using the AMBER99SB-ILDN force field. A total of 24 independent 500-ns production runs were performed.

In addition, potential of mean force (PMF) calculations were performed for **IL7Ra-55, IL7Ra-4H-01, and IL7Ra-3H-01 at pH 2** using the accelerated weight histogram (AWH) method.

## Systems and Simulation Conditions

| Protein variant | pH 2 | pH 7 | Replicates | Production length |
|----------------|------|------|------------|-------------------|
| IL7Ra-4H-01    | ✓    | ✓    | 3 | 500 ns each |
| IL7Ra-3H-01    | ✓    | –    | 3 | 500 ns each |
| IL7Rα-55       | ✓    | –    | 3 | 500 ns each |
| IL7Ra-4H-RD40  | ✓    | ✓    | 3 | 500 ns each |
| IL7Ra-4H-RD30  | ✓    | ✓    | 3 | 500 ns each |

**Protonation states:**
- **pH 2:** Asp and Glu were set to neutral forms (ASH and GLH, respectively), all His to doubly protonated state (HIP), C-terminus to COOH, and N-terminus retained as NH₃⁺.
- **pH 7:** Default AMBER99SB-ILDN protonation states were adopted.

All systems were solvated in a dodecahedral box with a minimum distance of 1.0 nm between protein and box edges, using TIP3P water. Na⁺ and Cl⁻ ions were added to a final concentration of 0.15 M.

### PMF Simulations (pH 2 only)

PMF calculations were performed for the following proteins at pH 2:
- IL7Ra-55 (3 target sites)
- IL7Ra-4H-01 (15 target sites)
- IL7Ra-3H-01 (9 target sites)

For each target site, **three independent 50-ns AWH runs** were performed.

PMF conditions:
- PMF profiles were computed using the Accelerated Weight Histogram (AWH) method as implemented in GROMACS 2024.6.
- The reaction coordinate was defined as the distance between the protein and the residue of interest.
- The coordinate was sampled over 0.0–2.0 nm using AWH with a force constant of 1000 kJ mol⁻¹ nm⁻².
- All PMF simulations were conducted in the NPT ensemble at 300 K and 1 bar, with temperature and pressure maintained by the V-rescale thermostat and the Parrinello–Rahman barostat, respectively.
- For each buried potential cleavage site, three independent 50-ns simulations were performed with randomly assigned initial velocities.
- The final PMF values were extracted from the well-converged last 1 ns (49–50 ns) of each replicate and are reported as the mean ± standard deviation (SD).
- The PMF was set to zero at 2.0 nm.

## Initial Structure Preparation

The initial structures of all variants were derived from AlphaFold2 (AF2) predicted models.

Because pH 2 and pH 7 require different protonation states, the AF2 models were processed separately for each pH condition **before** energy minimization:

- **pH 2:** The AF2 model was processed with manually assigned protonation states (Asp→ASH, Glu→GLH, His→HIP, C-terminus→COOH, N-terminus→NH₃⁺) using `pdb2gmx` (or equivalent), yielding a processed `.gro` file for each variant.
- **pH 7:** The AF2 model was processed using the default AMBER99SB-ILDN protonation states, yielding a separate processed `.gro` file for each variant.

These processed `.gro` files serve as the starting structures for all three replicates (rep1, rep2, rep3) under the corresponding pH condition. They are provided in the repository as:

01_Input_Files/<variant>/pH2/<variant>_pH2_processed.gro
01_Input_Files/<variant>/pH7/<variant>_pH7_processed.gro (if applicable)


The original AF2-predicted models (before protonation) are provided in `01_Input_Files/<variant>/common/` for reference.

## Software and Force Field

- **GROMACS:** 2024.6
- **Force field:** AMBER99SB-ILDN
- **Water model:** TIP3P
- **Trajectory analysis:** MDAnalysis 2.7.0, MDTraj, NumPy, pandas, Matplotlib
- **Python:**  3.10

## Simulation Protocol

- **Energy minimization:** Steepest descent, maximum 50,000 steps, convergence criterion 1000 kJ/mol/nm, position restraints on protein heavy atoms.
- **Equilibration:**
  - 200 ps NVT at 310 K using v-rescale thermostat.
  - 500 ps NPT at 310 K and 1 bar using C-rescale barostat.
  - Coordinates and energies recorded at 1-ps intervals.
- **Production:**
  - Three independent 500-ns runs per system with distinct initial random velocity seeds.
  - Temperature and pressure maintained at 310 K and 1 bar using v-rescale thermostat and Parrinello–Rahman barostat, respectively.
  - Coordinates saved every 10 ps (50,000 frames per replicate).
  - LINCS constraints on hydrogen bonds, 2-fs timestep.
  - Verlet cutoff scheme, 1.2 nm cutoff for Coulomb and van der Waals, PME for long-range electrostatics.

## Analyses and Scripts

All custom analysis scripts are located in `03_Analysis_Scripts/`. The following analyses were performed:

- **RMSD:** Backbone Cα RMSD relative to the corresponding AF2-predicted initial structure, computed using `gmx rms`. Time series were downsampled to every 50th frame and smoothed with a centered 50-frame moving average. (No separate custom script; directly using GROMACS built-in tools).
- **Dynamic cross-correlation (DCCM):** Cα-based DCCM using MDAnalysis. Only residue pairs separated by >10 residues were retained; |Cᵢⱼ| > 0.5 classified as significantly coupled. Script: `DCCM.py`.
- **Cleavage site solvent exposure:** Side-chain rSASA > 0.3 for pepsin, trypsin, and chymotrypsin sites, computed with an in-house Python script based on the Shrake–Rupley algorithm (MDTraj) and normalized using Miller side-chain max SASA values. Script: `sidechain_sasa_miller.py`.
- **Secondary structure (DSSP):** Using `gmx dssp`; trajectories downsampled every 20th frame. Script: `DSSP.py`.
- **PMF (AWH):** Accelerated Weight Histogram method in GROMACS 2024.6 for buried cleavage sites. Reaction coordinate: protein–residue distance 0.0–2.0 nm, force constant 1000 kJ mol⁻¹ nm⁻². Three independent 50-ns simulations per site. (No separate custom script; directly using GROMACS AWH module). Script: `PMF-draw.py`.Processed PMF data are provided in `04_Processed_Data/<variant>/pH2/pmf/`.
- **Core packing density:** Side-chain SASA < 10 Å² defines core residues; VvdW and Vtotal from `gmx sasa` and `gmx trjconv`. Script: `core-packing-density.sh`.
- **Interhelical hydrogen bonds:** `gmx hbond` with donor–acceptor distance cutoff 0.35 nm and D–H···A angle ≥ 150°. Helix pairs A–D analyzed. Script: `interhelical-hydrogen-bonds.sh`.

## Repository Structure

- `01_Input_Files/` – GROMACS input files for each protein variant and pH condition (topology, `.mdp`, processed starting structures, PMF inputs).
- `02_GROMACS_Key_Outputs/` – Key GROMACS outputs for each replicate (`.tpr`, final `.gro`, `md-log.7z`, `md_energy.csv`).
- `03_Analysis_Scripts/` – Python and shell scripts for all analyses (RMSD, DCCM, SASA, DSSP, PMF, packing density, hydrogen bonds).
- `04_Processed_Data/` – Processed data for figures (RMSD, DCCM, SASA, DSSP, packing density, hydrogen bonds, PMF).

Each protein variant and pH condition has its own subdirectory within these folders. The processed `.gro` files for each pH are shared by all three replicates. PMF data are only available for IL7Ra-55, IL7Ra-4H-01, and IL7Ra-3H-01 at pH 2.


**Note:** Each protein variant and pH condition has its own subdirectory. The processed `.gro` files for each pH are shared by all three replicates. PMF data are only available for IL7Ra-55, IL7Ra-4H-01, and IL7Ra-3H-01 at pH 2.

## Data Availability

All input files (`.mdp`, `.top`, `.itp`, processed `.gro`), `.tpr` files for exact reproduction, final structures (`.gro`), analysis scripts, and processed data are available in this GitHub repository: https://github.com/ShengquanLiu/TASCO/.

Raw energy files (`.edr`) and trajectories (`.xtc`) are too large to be shared. Key energy metrics extracted from the `.edr` files are provided as `md_energy.csv` in the repository. The complete `md.log` files have been compressed as `md-log.7z` and are included in the repository. The full `.edr` and `.xtc` files are available upon reasonable request.

## Citation

This repository accompanies a manuscript currently under review. Citation information will be updated upon publication. If you use these data or scripts before publication, please cite this GitHub repository:

Shengquan Liu et al. "De Novo Design of Protease-Resistant Minibinders via Topological and Sequential Co-Optimization." GitHub repository, 2026. https://github.com/ShengquanLiu/TASCO/

## License

This repository is provided for peer-review purposes. A formal license will be applied upon publication. For questions regarding reuse before publication, please contact the corresponding author.

## Contact

For questions, please contact:sqliu2016@sinano.ac.cn; kfpu2011@sinano.ac.cn; jli2006@sinano.ac.cn