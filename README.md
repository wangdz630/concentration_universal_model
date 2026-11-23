# Universal Models for Dark Matter Halos (Papers I, II, & III)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A unified, modular Python package linking **cosmological initial conditions**, **halo assembly histories (MAHs)**, and **halo internal density profiles** across cold dark matter (CDM), warm dark matter (WDM), and non-standard cosmologies.

This package integrates the complete three-part theoretical framework:
1. **Paper I (Population-Level Statistics)**: Universal fitting formulae for characteristic formation-time proxies $D(z_{\rm f})$ and statistical halo concentration relations (`peak`, `mean`, `median`) across CDM and $\Lambda$WDM cosmologies.
2. **Paper II (Individual Concentrations & Assembly Diversity)**: Predicting individual halo concentrations from assembly histories ($\omega_{\rm f}$) and modeling concentration distributions via simulation-calibrated Generalized Gamma distributions $P(\omega_{\rm f})$.
3. **Paper III (Universal MAHs & Bidirectional Archaeology)**: Parameter-free universal mass accretion histories $M(\tau)/M_0$ in self-similar $\tau$-space, forward prediction of continuous concentration evolution $c(z)$, and Bayesian Maximum A Posteriori (MAP) archaeological reconstruction of past MAHs from snapshot $(M_{\rm obs}, c_{\rm obs}, z_{\rm obs})$ observations.

Users can run each paper's model **completely independently**, or **couple them together** to explore structural evolution across cosmic time.

---

## Repository Structure

```text
.
├── README.md                           # Documentation
├── LICENSE                             # MIT License
├── .gitignore                          # Git ignore rules
│
├── cosmology_params.py                 # Cosmological parameter registry (Planck, WMAP, WDM, etc.)
├── GrowthFactor.py                     # Linear growth factor D(z) ODE solver
├── VarianceCalculator.py               # Mass variance sigma(M) via CAMB (top-hat & sharp-k)
├── halo_variance.py                    # Standalone variance & D(z_f) precomputation script
│
├── halo_model.py                       # Paper I & Paper II concentration models
├── mah_model.py                        # Paper III universal MAH, forward c(z), and Bayesian MAP inversion
├── halo_pipeline.py                    # Unified master CLI runner across all three papers
│
├── input_masses.dat                    # Sample input: halo masses [and optional z_f]
├── input_mass_concentration.dat        # Sample input: observed (M_obs, c_obs, z_obs) for reconstruction
├── input_merger_tree_sample.dat        # Sample input: progenitor trees (halo_id, z, M) for forward c(z)
│
├── output_halo_properties/             # Paper I outputs (statistical c: peak, mean, median)
├── output_individual_halo_properties/  # Paper II outputs (individual c & c(M) distributions)
└── output_mah_properties/             # Paper III outputs (MAH tracks, forward c(z), Bayesian MAP)
```

---

## Installation & Requirements

Ensure you have Python 3.8+ and the following scientific libraries installed:

```bash
pip install numpy scipy camb
```

Clone the repository:
```bash
git clone https://github.com/wangdz630/concentration_universal_model.git
cd concentration_universal_model
```

---

## Quick Start: Command Line Interface (CLI)

The package provides two convenient CLI entry points:
- **`halo_pipeline.py`**: Unified one-stop pipeline to run any paper or all papers in sequence.
- **`mah_model.py`**: Specialized, highly flexible tool for Paper III MAH modeling, forward tracking, and Bayesian inverse reconstruction.

---

### 1. Unified Pipeline (`halo_pipeline.py`)

Syntax:
```bash
python halo_pipeline.py [cosmology] [redshift] [mode] [halo_def]
```
- **`halo_def`** (optional): `vir`, `200c`, or `both` (default: runs **both `200c` and `vir`**).
- Outputs are saved into their respective clean directories.

#### Examples:

```bash
# -------------------------------------------------------------------------
# Paper I: Statistical concentration relations (saved to ./output_halo_properties/)
# -------------------------------------------------------------------------
python halo_pipeline.py planck18 0.0 peak            # Peak concentration (both 200c and vir)
python halo_pipeline.py planck18 0.0 peak vir        # Peak concentration (vir only)
python halo_pipeline.py planck18 0.0 mean            # Mean concentration
python halo_pipeline.py planck18 0.0 median          # Median concentration
python halo_pipeline.py planck18 0.0 all             # Run peak, mean, and median

# -------------------------------------------------------------------------
# Paper II: Individual concentrations & distributions (saved to ./output_individual_halo_properties/)
# -------------------------------------------------------------------------
python halo_pipeline.py planck18 0.0 relaxed         # NFW relaxed halos (both 200c & vir)
python halo_pipeline.py planck18 0.0 relaxed vir     # NFW relaxed halos (vir only)
python halo_pipeline.py planck18 0.0 all_halos       # NFW all halos (both 200c & vir)
python halo_pipeline.py planck18 0.0 relaxed_einasto # Einasto relaxed halos (200c)
python halo_pipeline.py planck18 0.0 individual      # Run all Paper II relations

# -------------------------------------------------------------------------
# Paper III: Universal MAH & Structural Evolution (saved to ./output_mah_properties/)
# -------------------------------------------------------------------------
python halo_pipeline.py planck18 0.0 mah             # Universal MAH tracks M(z)/M0 (both 200c & vir)
python halo_pipeline.py planck18 0.0 mah vir         # Universal MAH tracks (vir only)
python halo_pipeline.py planck18 0.0 c_evolution     # Continuous concentration evolution c(z) (relaxed & all)
python halo_pipeline.py planck18 0.0 c_evolution vir # Continuous c(z) (vir only)
python halo_pipeline.py planck18 0.0 reconstruct     # Bayesian MAP reconstruction from (M, c) (relaxed & all)
python halo_pipeline.py planck18 0.0 reconstruct vir # Bayesian MAP reconstruction (vir only)

# -------------------------------------------------------------------------
# Run All Papers Simultaneously
# -------------------------------------------------------------------------
python halo_pipeline.py planck18 0.0 all_papers      # Runs peak (I), relaxed (II), mah (III), and c_evolution (III)
```

---

### 2. Standalone Paper III CLI (`mah_model.py`)

`mah_model.py` provides dedicated subcommands for Paper III tasks with flexible inputs.

#### Input Data Options for MAHs and Concentrations:
- **Option A (Simulation Data Provided)**:
  - **Merger trees provided (`--tree file.dat`)**: The code measures the main-progenitor history and half-mass formation time $z_{\rm f}$ directly from the simulated tree.
  - **Formation redshifts provided (`Mass  z_f` in input file)**: The code uses the simulated $z_{\rm f}$ to evaluate the exact $\omega_{\rm f}$ and accretion history.
- **Option B (Simulation Data NOT Provided — Tree-Free Prediction)**:
  - **Only halo masses provided (`Mass` in input file)**: The code uses our **parameter-free Universal MAH model**, sampling the assembly parameter $\omega_{\rm f}$ from Paper II's simulation-calibrated distribution $P(\omega_{\rm f})$ (default `--mode sim`), theoretical EPS (`--mode eps`), or the population median (`--mode median`).

#### Task 1: Mass Accretion Histories (`mah`)
Compute individual or median MAH trajectories $M(z)/M_0$:
```bash
# Tree-free prediction using universal MAH model + Paper II P(omega_f):
python mah_model.py mah --input input_masses.dat --cosmology planck18 --halo-def vir

# If input file has 2 columns (Mass  z_f), it automatically uses the simulated z_f:
python mah_model.py mah --input input_masses.dat --mode sim_zf

# Driven directly by a simulated merger tree:
python mah_model.py mah --tree input_merger_tree_sample.dat --cosmology planck18

# Purely analytical population median MAHs (Paper III Figure 3):
python mah_model.py mah --input input_masses.dat --mode median
```

#### Task 2: Forward Concentration Evolution $c(z)$ (`c_evolution`)
Predict continuous concentration evolution along the progenitor branch from $z=5$ to $z=0$:
```bash
# Tree-free forward c(z) driven analytically by the universal MAH model:
python mah_model.py c_evolution --input input_masses.dat --cosmology planck18 --halo-def vir --relation relaxed

# Driven by an exact simulated merger tree (columns: halo_id  z  M(z)):
python mah_model.py c_evolution --tree input_merger_tree_sample.dat --cosmology planck18 --halo-def vir
```

#### Task 3: Bayesian Archaeological Inverse Reconstruction (`reconstruct`)
Reconstruct continuous historical MAHs and formation redshifts from snapshot observables $(M_{\rm obs}, c_{\rm obs}, [z_{\rm obs}])$:
```bash
# Reconstruct MAHs for observed clusters/halos:
python mah_model.py reconstruct --input input_mass_concentration.dat --cosmology planck18 --halo-def vir --relation relaxed
python mah_model.py reconstruct --input input_mass_concentration.dat --cosmology planck18 --halo-def 200c --relation all_halos
```

---

## Python API Usage

You can import `UniversalHaloModel` from `mah_model.py` in your own scripts or notebooks:

```python
import numpy as np
from mah_model import UniversalHaloModel

# Initialize model for Planck 2018 cosmology under virial definition
model = UniversalHaloModel('planck18', halo_def='vir')

M = np.array([1e12, 1e13, 1e14])  # Msun/h

# 1. Paper I: Statistical Concentration Relations
c_peak = model.concentration_paper1(M, z=0.0, statistic='peak')
print("Paper I Peak concentrations:", c_peak)

# 2. Paper II: Individual Halo Concentration
omega_f = np.array([1.10, 1.35, 1.70])  # Scaled formation parameter
c_indiv = model.concentration_paper2(M, omega_f, z=0.0, relation='relaxed')
print("Paper II Individual concentrations:", c_indiv)

# 3. Paper III: Universal Mass Accretion History (MAH)
z_eval = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
# Track M(z)/M0 for an individual halo with known half-mass redshift z_f = 1.2
ratio = model.mah(z_eval, M=1e14, z_f=1.2)
print("MAH M(z)/M0:", ratio)

# 4. Paper III: Forward Concentration Evolution c(z)
# Tree-free continuous c(z) trajectory driven by the universal MAH model:
c_z = model.concentration_track_model(z_eval=[0, 0.5, 1, 2, 3], M=1e14, omega_f=1.35)
print("Forward c(z):", c_z)

# 5. Paper III: Bayesian MAP Inverse Reconstruction
# Given an observed galaxy cluster: M_obs = 1e15 Msun/h, c_obs = 5.2 at z_obs = 0.44
recon = model.reconstruct(M_obs=1e15, c_obs=5.2, z_obs=0.44, z_eval=[0.5, 1.0, 2.0])
print("MAP inferred formation parameter omega_f:", recon['omega_f_MAP'])
print("MAP inferred formation redshift z_f:", recon['z_f_MAP'])
print("Historical MAH trajectory M(z)/M_obs:", recon['tracks'])
```

---

## Output Files & Column Reference

All outputs are saved as intuitive single text tables (one row per halo) with numbered column headers:

### 1. Paper I Outputs (`./output_halo_properties/`)
`halo_properties_<z>_<cosmo>_<def>_<mode>.dat`:
- `0:Mvir`: Halo mass $[h^{-1}M_\odot]$
- `1:sigma_M`: $\sigma(M, z=0)$
- `2:sigma_Mf`: $\sigma(0.5M, z=0)$
- `3:z`: Redshift
- `4:D_z`: Linear growth factor $D(z)$
- `5:D_zf_simfit`: Characteristic formation factor $D(z_{\rm f})$ from Paper I fit
- `8:c_simfit`: Predicted concentration (`peak`, `mean`, or `median`)
- `11:M_hm`: Half-mode mass $[h^{-1}M_\odot]$

### 2. Paper II Outputs (`./output_individual_halo_properties/`)
`c_distribution_<z>_<cosmo>_<def>_<mode>.dat`:
- `0:Mvir`: Halo mass $[h^{-1}M_\odot]$
- `5:omega_f_sim`: Scaled formation variable $\omega_{\rm f}$ (from input $z_{\rm f}$ or sampled from $P(\omega_{\rm f})$)
- `8:c_sim_eq5`: Exact multi-linear regression prediction
- `9:c_sim_mc`: Monte Carlo realization with intrinsic scatter $\mathcal{N}(0, \sigma_{\rm res}^2)$
- `10:c_eps_mc`: Realization using standard EPS $\omega_{\rm f}$

### 3. Paper III Outputs (`./output_mah_properties/`)
- **Task 1: MAH Predictions** (`mah_<z>_<cosmo>_<def>_<mode>.dat`):
  - `0:M0[Msun/h]`: Descendant / observed halo mass
  - `1:z_f`: Half-mass formation redshift $z_{\rm f}$
  - `2:omega_f`: Scaled formation variable $\omega_{\rm f}$
  - `3:beta`: MAH shape parameter $\beta = A_{\rm hm}\omega_{\rm f} - \omega_{\rm f,peak}$
  - `4..N:M(z)/M0`: Normalized mass growth ratios $M(z)/M_0$ at evaluated redshifts

- **Task 2: Forward Concentration Evolution** (`c_evolution_<z>_<cosmo>_<def>_<relation>_<src>.dat`):
  - `0:M0[Msun/h]`: Descendant halo mass
  - `1:z_f`: Half-mass formation redshift
  - `2..N:c(z)`: Continuous concentration evolution $c(z)$ along the main progenitor branch

- **Task 3: Bayesian Reconstructed MAHs** (`reconstruct_<cosmo>_<def>_<relation>.dat`):
  - `0:M_obs[Msun/h]`, `1:c_obs`, `2:z_obs`: Input observables
  - `3:z_f_MAP`: Bayesian MAP inferred half-mass formation redshift
  - `4:omega_f_MAP`: Bayesian MAP inferred formation parameter
  - `5..N:M(z)/M_obs`: Reconstructed past MAH growth trajectory into the past

---

## Supported Cosmological Models

Configurations defined in `cosmology_params.py`:
- `planck18`, `planck15`, `planck13`
- `WMAP1`, `WMAP3`, `WMAP5`, `WMAP7`, `WMAP9`
- `illustris`, `bolshoi`, `multidark-planck`, `millennium`
- `SCDM` (Einstein--de Sitter), `OCDM` (Open CDM), `wCDM` (Dynamical Dark Energy)
- `LWDM` (Warm Dark Matter with thermal relic masses, e.g. $m_{\rm WDM} = 1.5\,\mathrm{keV}$)
- `scale-free` (Power-law spectra with $n = -2$)

---

## Citations

If you use this code package in your research, please cite our papers:

1. **Paper I (Peak & Statistical Concentration)**:  
   Wang, D.-Z., Lin, W., & Luan, T.-C. 2026, *The Astrophysical Journal*, 1006, 61.  
   DOI: [10.3847/1538-4357/ae81a2](https://doi.org/10.3847/1538-4357/ae81a2)

2. **Paper II (Individual Halo Concentration & Distributions)**:  
   Wang, D.-Z., Lin, W., Zhu, W., & Luan, T.-C. 2026, *The Astrophysical Journal*, submitted.  
   arXiv: [2611.17676](https://arxiv.org/abs/2611.17676)

3. **Paper III (Universal MAH & Structural Evolution)**:  
   Wang, D.-Z., Lin, W., Zhu, W., & Luan, T.-C. 2026, *The Astrophysical Journal*, submitted.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
