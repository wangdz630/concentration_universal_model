# Halo Concentration Model

A Python package for computing halo mass variance, formation-time proxies, and halo concentration predictions using a universal model calibrated across multiple cosmologies, halo definitions, and concentration statistics.

This repository implements a two-step workflow:

1. compute halo mass variance from an input halo mass list;
2. predict halo formation-time quantities and halo concentration using the universal model.

The package supports both standard CDM-like cosmologies and LWDM cosmologies with an analytic half-mode mass correction.

---

## Features

- Compute halo variance from input halo masses
- Predict halo formation-time proxy `D(z_f)`
- Predict halo concentration using a universal calibrated model
- Support for multiple cosmologies:
  - `ΛCDM`
  - `OCDM`
  - `wCDM`
  - `SCDM`
  - `LWDM`
  - scale-free / EdS-like setups
- A wide range of cosmological models is provided in `cosmology_params.py`; users may also define additional cosmologies of interest by extending this file.
- Support for multiple halo definitions:
  - `vir`
  - `200c`
- Support for multiple concentration relations:
  - `peak`
  - `mean`
  - `median`
- Automatic LWDM half-mode mass calculation from cosmological parameters

---

## Repository Structure

```text
.
├── README.md
├── input_masses.dat
├── cosmology_params.py
├── GrowthFactor.py
├── VarianceCalculator.py
├── halo_variance.py
├── halo_model.py
├── halo_pipeline.py
└── output_halo_properties/
```

## Main Files

- `cosmology_params.py`  
  Cosmological parameter configurations.

- `GrowthFactor.py`  
  Linear growth factor calculation.

- `VarianceCalculator.py`  
  Halo variance calculation tools.

- `halo_variance.py`  
  Script to compute halo variance from an input mass file.

- `halo_model.py`  
  Universal halo formation and concentration model.

- `halo_pipeline.py`  
  End-to-end pipeline for variance and concentration prediction.

---

## Quick Start

Prepare an input file named:

```text
input_masses.dat
```

Then run the full pipeline, for example:

```bash
python halo_pipeline.py planck18 0.0 
python halo_pipeline.py planck18 0.0 median
```

If you only want to compute the halo mass variance, for example:

```bash
python halo_variance.py input_masses.dat planck18 0.0
```

---

## Model Description

The universal formation-time model is

```math
D(z_f)_{\rm CDM}
=
D(z)\left[
a_1 \left(\frac{M}{10^{12}}\right)^{b_1}
+
a_2 D(z)^{b_2}
\right].
```

For LWDM cosmologies, a half-mode mass correction is applied:

```math
D(z_f)
=
D(z_f)_{\rm CDM}
\left[
1 + \eta \left(\frac{M_{\rm hm}}{M}\right)^\mu
\right].
```

The concentration model is

```math
c = \alpha \nu_{\rm eff}^{\beta} + \gamma,
```

with

```math
\nu_{\rm eff} = \nu_{\rm f}\left(1 + A \frac{M_{\rm hm}}{M}\right),
```

and

```math
\nu_{\rm f}
=
\frac{\delta_c}{\sigma(M)\,D(z)\,D(z_f)}.
```

The implemented model supports different halo definitions (`vir`, `200c`) and different relations (`peak`, `mean`, `median`).

---

## Scientific Scope

This package is intended for:

- halo concentration modeling
- halo formation-time proxy studies
- cross-cosmology comparisons
- CDM and LWDM halo structure analysis

It is designed as a research implementation of the universal concentration model used in this work.

---

## Citation

If you use this code in scientific work, please cite the associated paper of this repository.
