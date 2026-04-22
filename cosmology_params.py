"""
Cosmological Parameter Configurations
=====================================
This module provides pre-defined cosmological parameter sets for various
cosmological models, including ΛCDM, wCDM, SCDM, OCDM, scale-free, and WDM.
All parameters are defined at the present epoch (z=0) unless otherwise specified.

Below cosmology params from paper: COLOSSUS: A Python Toolkit for Cosmology, Large-scale Structure, 
and Dark Matter Halos (Benedikt Diemer)
"""

cosmology_configs = {
    'planck18': {
        'h': 0.6766,
        'omega_m': 0.3111,
        'omega_b': 0.0490,
        'n': 0.9665,
        'sigma_8': 0.8102,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'planck15': {
        'h': 0.6774,
        'omega_m': 0.3089,
        'omega_b': 0.0486,
        'n': 0.9667,
        'sigma_8': 0.8159,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'planck13': {
        'h': 0.6777,
        'omega_m': 0.3071,
        'omega_b': 0.0483,
        'n': 0.9611,
        'sigma_8': 0.8288,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'WMAP9': {
        'h': 0.6932,
        'omega_m': 0.2865,
        'omega_b': 0.0463,
        'n': 0.9608,
        'sigma_8': 0.8200,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'WMAP7': {
        'h': 0.7020,
        'omega_m': 0.2743,
        'omega_b': 0.0458,
        'n': 0.9680,
        'sigma_8': 0.8160,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'WMAP5': {
        'h': 0.7050,
        'omega_m': 0.2732,
        'omega_b': 0.0456,
        'n': 0.9600,
        'sigma_8': 0.8120,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'WMAP3': {
        'h': 0.7350,
        'omega_m': 0.2342,
        'omega_b': 0.0413,
        'n': 0.9510,
        'sigma_8': 0.7420,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'WMAP1': {
        'h': 0.7200,
        'omega_m': 0.2700,
        'omega_b': 0.0463,
        'n': 0.9900,
        'sigma_8': 0.9000,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'illustris': {
        'h': 0.7040,
        'omega_m': 0.2726,
        'omega_b': 0.0456,
        'n': 0.9630,
        'sigma_8': 0.8090,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'bolshoi': {
        'h': 0.7000,
        'omega_m': 0.2700,
        'omega_b': 0.0469,
        'n': 0.9500,
        'sigma_8': 0.8200,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'multidark-planck': {
        'h': 0.6780,
        'omega_m': 0.3070,
        'omega_b': 0.0480,
        'n': 0.9600,
        'sigma_8': 0.8290,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'millennium': {
        'h': 0.7300,
        'omega_m': 0.2500,
        'omega_b': 0.0450,
        'n': 1.0000,
        'sigma_8': 0.9000,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    }, 
    'SCDM': {
        'h': 0.50,
        'omega_m': 1.0,
        'omega_b': 0.048,
        'n': 1.0,
        'sigma_8': 0.55,
        'delta_c': 1.686,
        'cosmology_type': 'CDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
    'OCDM': {
        'h': 0.667,
        'omega_m': 0.30,
        'omega_b': 0.048,
        'n': 1.0,
        'sigma_8': 1.0,
        'delta_c': 1.686,
        'cosmology_type': 'OCDM',
        'w_0': -1.0,
        'w_a': 0.0
    },
## wCDM by using planck18 parameter
    'wCDM': {
        'h': 0.6766,
        'omega_m': 0.3111,
        'omega_b': 0.0490,
        'n': 0.9665,
        'sigma_8': 0.8102,
        'delta_c': 1.686,
        'cosmology_type': 'wCDM',
        'w_0': -0.7,
        'w_a': 0.0
    },
## LWDM by using planck18 parameter, m_kev=1.5 keV
    'LWDM': {
        'h': 0.6766,
        'omega_m': 0.3111,
        'omega_b': 0.0490,
        'n': 0.9665,
        'sigma_8': 0.8102,
        'delta_c': 1.686,
        'cosmology_type': 'LWDM',
        'w_0': -1.0,
        'w_a': 0.0,
        'm_keV': 1.5
    },
## SF simulations with spectral index n
    'scale-free': {
        'h': 0.6766,
        'omega_m': 1.0,
        'omega_b': 0.0,
        'n': 1.0,
        'sigma_8': 1.0,
        'delta_c': 1.686,
        'cosmology_type': 'EdS',
        'w_0': -1.0,
        'w_a': 0.0,
    }
}
