# mace_phonons

This repository contains phonon calculations and evaluations for the MACE MP-0b3 model.

### dft_pbe_phonons

The phonon calculations for the benchmark set were performed using density functional theory (DFT) with the PBE functional. These calculations serve as reference data for evaluating the performance of the MACE model.  
The calculations were carried out using the finite-displacement method with `Phonopy` ([link](https://phonopy.github.io/phonopy/)) and the Atomate2 `PhononMaker` ([link](https://github.com/materialsproject/atomate2)).  
For the benchmark set, mainly chalcogenide-based and thermoelectric materials were selected. In the directory `dft_pbe_phonons/{mp_id}`, you can find the `phonopy.yaml` files for the different materials, with the `mp-id` as used in the Materials Project ([link](https://next-gen.materialsproject.org)).

#### DFT Settings:
The DFT calculations were performed with VASP.    
The geometric optimization was performed twice with the following settings:     
`ALGO` = Normal    
`EDIFF` = 1e-07      
`EDIFFG` = 1e-06     
`ENCUT` = 520    
`IBRION` = 2    
`ISIF` = 3  
`ISMEAR` = 0    
`ISPIN` = 1     
`LASPH` = True      
`LORBIT` = 11   
`LREAL` = False     
`NELM` = 100      
`NPAR` = 8  
`NSW` = 99  
`PREC` = Accurate   
`SIGMA` = 0.05  
The number of k-points was set to 300 $\frac{\text{k-points}}{Å^{-3}}$.

The single atom displaced supercells were constructed with a minimum edge length of 20 Å and a cubic shape of 90° angles.
Single atoms were displaced by 0.01 Å. 
The ionic convergence criterion was set to 1e−6 eV and the calculation was performed using only a single k-point at $\Gamma$.  

### evaluation

The script used to compare the DFT-PBE phonons and MP-0b3 phonons is located at `evaluation/calculate_phonon_properties.py`.  
The results can be found in the file `phonon_properties.csv`.

### functions

Plotting and metric functions used during evaluation.

### mp_0b3_medium_phonons

Phonon calculations with the MP-0b3 model, performed using the Atomate2 `PhononMaker`.
Geometric optimisation of the primitive cell and force calculation of the single atom displaced supercells were performed with the MP-0b3 model.
During relaxation the stop criterion was set to 1e-5 eV/Å.
The same supercell matrix was used for the DFT and MACE calculations.

### plots
DFT phonon band structures (black) and MACE-MP-0b3 phonon band structures (red) plotted together.

### phonon_properties.csv

The columns in the `phonon_properties.csv` file have the following meanings:

- `rmse`: Root Mean Square Error between DFT and MACE phonons along the selected k-path.
- `mae`: Mean Absolute Error between DFT and MACE phonons along the selected k-path.
- `r2`: R² score between DFT and MACE phonons along the selected k-path.
- `phonon_rmse`: Root Mean Square Error between DFT and MACE phonons evaluated over the entire Brillouin zone.
- `phonon_mae`: Mean Absolute Error between DFT and MACE phonons over the entire Brillouin zone.
- `phonon_scaled_rmse`: DFT and MACE phonons are scaled to the same maximum frequency, then RMSE is computed over the Brillouin zone. This metric compares the shape of the phonon bands.
- `phonon_rrmse`: Relative RMSE, i.e., the RMSE between DFT and MACE phonons over the Brillouin zone normalized by the RMS of the DFT phonons.
- `rmsd`: Root Mean Square Distance between the DFT-relaxed primitive cell and the MACE-MP-0b3 relaxed primitive cell, as defined in `pymatgen.analysis.structure_matcher.StructureMatcher`.
- `maxd`: Maximum displacement between the DFT and MACE-MP-0b3 relaxed primitive cells, also from `StructureMatcher`.
- `spacegroup`: Space group of the MACE-MP-0b3 relaxed primitive cell.
- `min_frequency`: Minimum frequency of the MACE-MP-0b3 phonons.
- `max_frequency`: Maximum frequency of the MACE-MP-0b3 phonons.
- `mean_frequency`: Mean frequency of the MACE-MP-0b3 phonons.
- `volume`: Volume of the MACE-MP-0b3 relaxed primitive cell.
- `benchmark_spacegroup`: Space group of the DFT-PBE relaxed primitive cell.
- `benchmark_min_frequency`: Minimum frequency of the DFT-PBE phonons.
- `benchmark_max_frequency`: Maximum frequency of the DFT-PBE phonons.
- `benchmark_mean_frequency`: Mean frequency of the DFT-PBE phonons.
- `benchmark_volume`: Volume of the DFT-PBE relaxed primitive cell.

---

## A foundation model for atomistic materials chemistry

The MP-0b3 model is based on the [MACE foundation model](https://arxiv.org/abs/2401.00096) and will be included in a future update of that work together with the phonon data.


```bibtex
@article{batatia2023foundation,
  title={A foundation model for atomistic materials chemistry},
  author={Batatia, Ilyes and Benner, Philipp and Chiang, Yuan and Elena, Alin M and Kov{\'a}cs, D{\'a}vid P and Riebesell, Janosh and Advincula, Xavier R and Asta, Mark and Avaylon, Matthew and Baldwin, William J and others},
  journal={arXiv preprint arXiv:2401.00096},
  year={2023}
}