# Universal Models for Dark Matter Halos (Papers I, II, & III)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A unified, modular Python package linking **cosmological initial conditions**, **halo assembly histories (MAHs)**, and **halo internal density profiles** across cold dark matter (CDM), warm dark matter (WDM), and non-standard cosmologies.

This package integrates the complete three-part theoretical framework:
1. **Paper I (Population-level Statistics)**: Characteristic formation-time proxies and statistical halo concentration relations (`peak`, `mean`, `median`) across CDM and $\Lambda$WDM cosmologies.
2. **Paper II (Individual Halo Concentrations & Assembly Diversity)**: Predicting individual halo concentrations from assembly histories ($\omega_{\rm f}$) and modeling concentration distributions via simulation-calibrated Generalized Gamma distributions $P(\omega_{\rm f})$.
3. **Paper III (Universal Mass Accretion Histories & Bidirectional Archaeology)**: Parameter-free universal mass accretion histories $M(\tau)/M_0$ in self-similar $\tau$-space, forward prediction of continuous concentration evolution $c(z)$, and Bayesian Maximum A Posteriori (MAP) archaeological reconstruction of past MAHs from snapshot $(M_{\rm obs}, c_{\rm obs}, z_{\rm obs})$ observations.

Users can run each paper's model **completely independently**, or **couple them together** to explore structural evolution across cosmic time.

---

## Repository Structure

```text
.
├── README.md                           # Comprehensive documentation
├── LICENSE                             # MIT License
├── cosmology_params.py                 # Cosmological parameter registry (Planck, WMAP, WDM, etc.)
├── GrowthFactor.py                     # Linear growth factor D(z) differential equation solver
├── VarianceCalculator.py               # Mass variance sigma(M) via CAMB (top-hat & sharp-k)
├── halo_variance.py                    # Standalone variance & D(z_f) precomputation script
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

Ensure you have Python 3.8+ and the following dependencies installed:

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
- **`halo_pipeline.py`**: A one-stop pipeline for running any paper or all papers in sequence.
- **`mah_model.py`**: A specialized, highly flexible tool for Paper III MAH modeling, forward tracking, and Bayesian inverse reconstruction.

### 1. Unified Pipeline (`halo_pipeline.py`)

Syntax: `python halo_pipeline.py [cosmology] [redshift] [mode] [halo_def]`
- If `halo_def` is omitted, the pipeline runs **both `200c` and `vir`** by default.
- Specify `vir` or `200c` to run for a specific halo definition only.

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
python halo_pipeline.py planck18 0.0 c_evolution     # Continuous concentration evolution c(z)
python halo_pipeline.py planck18 0.0 c_evolution vir # Continuous c(z) (vir only)
python halo_pipeline.py planck18 0.0 reconstruct     # Bayesian MAP reconstruction from (M, c)
python halo_pipeline.py planck18 0.0 reconstruct vir # Bayesian MAP reconstruction (vir only)

# -------------------------------------------------------------------------
# Run All Papers Simultaneously
# -------------------------------------------------------------------------
python halo_pipeline.py planck18 0.0 all_papers      # Runs peak, relaxed, mah, c_evolution
python halo_pipeline.py planck18 0.0 all_papers vir  # All papers for vir only
```

---

### 2. Standalone Paper III CLI (`mah_model.py`)

`mah_model.py` provides dedicated subcommands for advanced Paper III tasks:

#### Task A: Mass Accretion Histories (`mah`)
Compute individual or median MAH trajectories $M(z)/M_0$:
```bash
# Automatically detects if input_masses.dat has 2 columns (Mass z_f) or 1 column (Mass)
python mah_model.py mah --input input_masses.dat --cosmology planck18 --halo-def vir

# Generate purely analytical population median MAHs (Paper III Figure 3)
python mah_model.py mah --input input_masses.dat --omega median

# Generate synthetic Monte Carlo MAH ensembles drawn from Paper II's P(omega_f) (Paper III Figure 4)
python mah_model.py mah --input input_masses.dat --omega sim
```

#### Task B: Forward Concentration Evolution $c(z)$ (`c_evolution`)
Predict continuous concentration evolution from $z=5$ to $z=0$ (Paper III Section 4, Figure 7):
```bash
# Tree-independent forward c(z) driven by the universal MAH model:
python mah_model.py c_evolution --input input_masses.dat --cosmology planck18

# Driven by an exact simulated merger tree (columns: halo_id  z  M(z)):
python mah_model.py c_evolution --tree input_merger_tree_sample.dat --cosmology planck18
```

#### Task C: Bayesian Archaeological Inverse Reconstruction (`reconstruct`)
Reconstruct continuous historical MAHs and $1\sigma$ credible bands from snapshot observables $(M_{\rm obs}, c_{\rm obs}, z_{\rm obs})$ (Paper III Section 5, Figure 8 & Figure 9):
```bash
python mah_model.py reconstruct --input input_mass_concentration.dat --cosmology planck18 --halo-def vir
```

---

## Python API: Unified Python Usage

You can import `UniversalHaloModel` from `mah_model.py` to seamlessly execute calculations from Papers I, II, and III in your own scripts or Jupyter notebooks:

```python
import numpy as np
from mah_model import UniversalHaloModel

# Initialize model for Planck 2018 cosmology under virial definition
model = UniversalHaloModel('planck18', halo_def='vir')

# =========================================================================
# 1. Paper I: Statistical Concentration Relations
# =========================================================================
M = np.array([1e12, 1e13, 1e14])  # Msun/h
c_peak = model.concentration_paper1(M, z=0.0, statistic='peak')
print("Paper I Peak concentrations:", c_peak)

# =========================================================================
# 2. Paper II: Individual Halo Concentration
# =========================================================================
omega_f = np.array([1.10, 1.35, 1.70])  # Assembly parameter
c_indiv = model.concentration_paper2(M, omega_f, z=0.0, relation='relaxed')
print("Paper II Individual concentrations:", c_indiv)

# =========================================================================
# 3. Paper III: Universal Mass Accretion History (MAH)
# =========================================================================
z_eval = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
# Track M(z)/M0 for an individual halo with known half-mass redshift z_f = 1.2
ratio = model.mah(z_eval, M=1e14, z_f=1.2)
print("MAH M(z)/M0:", ratio)

# =========================================================================
# 4. Paper III: Forward Concentration Evolution c(z)
# =========================================================================
# Continuous c(z) trajectory (tree-free, driven by universal MAH model):
c_z = model.concentration_track_model(z_eval=[0, 0.5, 1, 2, 3], M=1e14, omega_f=1.35)
print("Forward c(z):", c_z)

# =========================================================================
# 5. Paper III: Bayesian MAP Inverse Reconstruction
# =========================================================================
# Given an observed galaxy cluster: M_obs = 1e15 Msun/h, c_obs = 5.2 at z_obs = 0.44
recon = model.reconstruct(M_obs=1e15, c_obs=5.2, z_obs=0.44, z_eval=[0.5, 1.0, 2.0])
print("MAP inferred formation parameter omega_f:", recon['omega_f_MAP'])
print("MAP inferred formation redshift z_f:", recon['z_f_MAP'])
print("Historical MAH trajectory M(z)/M_obs:", recon['tracks'])
```

---

## Theoretical Framework

### 1. Paper I — Statistical Concentration Model
Relates peak height $\nu$ and characteristic formation growth factor $D(z_{\rm f})$:
```math
D(z_{\rm f})_{\rm CDM} = D(z)\left[a_1 \left(\frac{M}{10^{12}\,h^{-1}M_\odot}\right)^{b_1} + a_2 D(z)^{b_2}\right],
```
```math
c = \alpha \nu_{\rm eff}^\beta + \gamma, \qquad \nu_{\rm eff} = \frac{\delta_{\rm sc}(0)}{\sigma(M) D(z) D(z_{\rm f})} \left(1 + A \frac{M_{\rm hm}}{M}\right).
```

### 2. Paper II — Individual Halo Concentration Relation
Decouples cosmic expansion, mass variance, individual assembly history $\omega_{\rm f}$, and WDM suppression:
```math
\log_{10} c_{\rm vir}(z) = A + \alpha_D \log_{10} D(z) + \alpha_\sigma \log_{10} \sigma(M) + \beta \log_{10}(1 + \omega_{\rm f}) + \gamma \log_{10} A_{\rm hm},
```
where the scaled formation variable is
```math
\omega_{\rm f} = \frac{\delta_{\rm sc}(z_{\rm f}) - \delta_{\rm sc}(z)}{\sqrt{\sigma^2(0.5M) - \sigma^2(M)}}.
```

### 3. Paper III — Universal Parameter-Free MAH Model
In dimensionless self-similar time $\tau(z) = [D(z)^{-1} - 1] / [D(z_{\rm f})^{-1} - 1]$, mass variance cancels identically:
```math
\frac{M(\tau)}{M_0} = (1 + \tau)^\beta 2^{-(1 + \beta)\tau}, \qquad \beta = A_{\rm hm} \omega_{\rm f} - \omega_{\rm f,peak},
```
where $\omega_{\rm f,peak} = 1.25$ (vir) and $1.29$ (200c), and $A_{\rm hm} = 1 + a (M_{\rm hm}/M_0)$ ($A_{\rm hm} = 1$ for CDM). Contains **zero free fitting parameters**.

### 4. Paper III — Bayesian MAP Archaeological Inversion
Inverting the linear concentration observable $y = \beta u + \epsilon$ (where $u \equiv \log_{10}(1 + \omega_{\rm f})$) with a Gaussian cosmological prior $u \sim \mathcal{N}(0.368, 0.094^2)$ yields the exact conjugate posterior:
```math
\mu_{\rm post} = \frac{w_{\rm prior}\mu_{\rm prior} + w_{\rm like}(y/\beta)}{w_{\rm prior} + w_{\rm like}}, \qquad \sigma_{\rm post} = \frac{1}{\sqrt{w_{\rm prior} + w_{\rm like}}}.
```

---

## Output Files & Column Reference

All outputs are saved as single, intuitive text tables (one row per halo / cluster) with explicit column number headers matching the Paper I and Paper II formatting convention.

### 1. Paper I Outputs (`./output_halo_properties/`)
- `0:Mvir`: Halo mass $[h^{-1}M_\odot]$
- `1:sigma_M`: $\sigma(M, z=0)$
- `2:sigma_Mf`: $\sigma(0.5M, z=0)$
- `3:z`: Redshift
- `4:D_z`: Linear growth factor $D(z)$
- `5:D_zf_simfit`: Characteristic formation factor $D(z_{\rm f})$ from Paper I fit
- `8:c_simfit`: Predicted concentration (`peak`, `mean`, or `median`)
- `11:M_hm`: Half-mode mass $[h^{-1}M_\odot]$

### 2. Paper II Outputs (`./output_individual_halo_properties/`)
- `0:Mvir`: Halo mass $[h^{-1}M_\odot]$
- `5:omega_f_sim`: Scaled formation variable $\omega_{\rm f}$ (from input $z_{\rm f}$ or sampled from $P(\omega_{\rm f})$)
- `8:c_sim_eq5`: Exact multi-linear regression prediction
- `9:c_sim_mc`: Monte Carlo realization with intrinsic scatter $\mathcal{N}(0, \sigma_{\rm res}^2)$
- `10:c_eps_mc`: Realization using standard EPS $\omega_{\rm f}$

### 3. Paper III Outputs (`./output_mah_properties/`)
- **Task 1: MAH Predictions** (`mah_predictions_<z>_<cosmo>_<def>_<mode>.dat`):
  - `0:M0[Msun/h]`: Present/observed halo mass
  - `1:z_obs`: Observation redshift
  - `2:omega_f`: Scaled formation variable $\omega_{\rm f}$
  - `3:z_f`: Half-mass formation redshift $z_{\rm f}$
  - `4:beta`: MAH shape parameter $\beta = A_{\rm hm}\omega_{\rm f} - \omega_{\rm f,peak}$
  - `5..N:M(z)/M0`: Normalized mass accretion ratios $M(z)/M_0$ at evaluated redshifts (e.g. $z=0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0$)

- **Task 2: Forward Concentration Evolution** (`c_evolution_<cosmo>_<def>_<relation>_<src>.dat`):
  - `0:M0[Msun/h]`: Descendant halo mass
  - `1..N:c(z)`: Continuous concentration evolution $c(z)$ evaluated from $z=0$ to high redshift along the progenitor branch

- **Task 3: Bayesian Reconstructed MAHs** (`reconstructed_mah_<cosmo>_<def>_<relation>.dat`):
  - `0:M_obs[Msun/h]`, `1:c_obs`, `2:z_obs`: Observed input parameters
  - `3:omega_f_MAP`: Bayesian MAP estimated formation parameter
  - `4:z_f_MAP`: Bayesian MAP inferred half-mass formation redshift
  - `5:M0_MAP[Msun/h]`: Descendant $z=0$ mass
  - `6..N:M(z)/M_obs`: Reconstructed past MAH growth trajectory $M(z)/M_{\rm obs}$ into the past

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

If you use this code package in your research, please cite the corresponding papers:

```bibtex
@ARTICLE{2026ApJ..1006...61W,
       author = {{Wang}, Dao-Zhou and {Lin}, Weipeng and {Luan}, Tian-Cheng},
        title = "{Universal Fitting Formulae for the Peak Concentration of Dark Matter Halos}",
      journal = {\apj},
     keywords = {Cosmology, Dark matter, Dark matter halos, N-body simulations},
         year = 2026,
       volume = {1006},
       number = {1},
          eid = {61},
        pages = {61},
          doi = {10.3847/1538-4357/ae81a2}
}

@ARTICLE{2026arXiv261117676W,
       author = {{Wang}, Dao-Zhou and {Lin}, Weipeng and {Zhu}, Weishan and {Luan}, Tian-Cheng},
        title = "{Predicting Individual Dark Matter Halo Concentrations from Their Assembly Histories}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2026,
        month = nov,
          eid = {arXiv:2611.17676},
        pages = {arXiv:2611.17676},
archivePrefix = {arXiv},
       eprint = {2611.17676}
}

@ARTICLE{PaperIII,
       author = {{Wang}, Dao-Zhou and {Lin}, Weipeng and {Zhu}, Weishan and {Luan}, Tian-Cheng},
        title = "{A Universal Model for the Mass Accretion Histories of Dark Matter Halos}",
      journal = {AAS Journals, submitted},
         year = 2026
}
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
