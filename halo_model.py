"""
Calculate halo concentrations using the new universal model:

D(z_f)_CDM = D(z) * [a1 * (M/1e12)^(b1) + a2 * D(z)^(b2)]
D(z_f) = D(z_f)_CDM * [1 + eta * (M_hm/M)^mu]

c = alpha * nu_eff^(beta) + gamma
nu_eff = nu_peak * (1 + A * M_hm/M)
nu_peak = delta_c / [sigma(M) * D(z) * D(z_f)]

This version handles both CDM and LWDM cosmologies, supporting 
multiple halo definitions (200c, vir) and relations (peak, mean, median).
"""

import numpy as np
import sys
import os

# =========================================================
# Import cosmology config
# =========================================================
try:
    from cosmology_params import cosmology_configs
except ImportError:
    print("Warning: cosmology_params module not found.")
    cosmology_configs = {}

# =========================================================
# Model Parameters Dictionary
# =========================================================
MODEL_PARAMS = {
    '200c': {
        'peak':   {'a1': 0.17, 'b1': 0.14, 'a2': 0.37, 'b2': -0.49, 'eta': 0.14, 'mu': 1.06, 'alpha': 8.74,  'beta': -1.04, 'gamma': 2.53, 'A': 3.33},
        'mean':   {'a1': 0.20, 'b1': 0.12, 'a2': 0.35, 'b2': -0.48, 'eta': 0.13, 'mu': 0.91, 'alpha': 9.47,  'beta': -1.22, 'gamma': 3.16, 'A': 1.78},
        'median': {'a1': 0.18, 'b1': 0.13, 'a2': 0.35, 'b2': -0.48, 'eta': 0.13, 'mu': 1.31, 'alpha': 9.18,  'beta': -1.14, 'gamma': 2.93, 'A': 2.59}
    },
    'vir': {
        'peak':   {'a1': 0.21, 'b1': 0.12, 'a2': 0.36, 'b2': -0.45, 'eta': 0.13, 'mu': 1.08, 'alpha': 11.68, 'beta': -1.17, 'gamma': 2.68, 'A': 3.92},
        'mean':   {'a1': 0.24, 'b1': 0.11, 'a2': 0.34, 'b2': -0.45, 'eta': 0.12, 'mu': 0.94, 'alpha': 12.26, 'beta': -1.26, 'gamma': 3.23, 'A': 2.26},
        'median': {'a1': 0.24, 'b1': 0.10, 'a2': 0.34, 'b2': -0.45, 'eta': 0.13, 'mu': 1.28, 'alpha': 11.97, 'beta': -1.19, 'gamma': 2.98, 'A': 3.03}
    }
}

# =========================================================
# Physics Functions
# =========================================================
def universal_D_zf(D_z, Mvir, params, M_hm=0.0):
    D_z = np.asarray(D_z, dtype=float)
    Mvir = np.asarray(Mvir, dtype=float)
    M12 = Mvir / 1.0e12
    D_zf_cdm = D_z * (params['a1'] * M12**params['b1'] + params['a2'] * D_z**params['b2'])

    if np.any(np.asarray(M_hm) > 0):
        return D_zf_cdm * (1.0 + params['eta'] * (M_hm / Mvir)**params['mu'])
    return D_zf_cdm

def D_zf_eps_func(D_z, sigma_M, sigma_M_f, omega_f=0.75, delta_c=1.686):
    D_z = np.asarray(D_z, dtype=float)
    sigma_M = np.asarray(sigma_M, dtype=float)
    sigma_M_f = np.asarray(sigma_M_f, dtype=float)
    Delta_S2 = np.maximum(sigma_M_f**2 - sigma_M**2, 0.0)
    Delta_S = np.sqrt(Delta_S2)
    return D_z / (1.0 + (omega_f / delta_c) * D_z * Delta_S)

def universal_concentration(D_z, D_zf, sigma_M, Mvir, params, M_hm=0.0, delta_c=1.686):
    D_z = np.asarray(D_z, dtype=float)
    D_zf = np.asarray(D_zf, dtype=float)
    sigma_M = np.asarray(sigma_M, dtype=float)
    Mvir = np.asarray(Mvir, dtype=float)

    denom = sigma_M * D_z * D_zf
    nu_peak = delta_c / denom
    nu_eff = nu_peak * (1.0 + params['A'] * M_hm / Mvir)
    return params['alpha'] * nu_eff**params['beta'] + params['gamma']

def calculate_lwdm_half_mode_mass(cosmo_params, m_wdm_keV=1.5):
    h = cosmo_params.get('h', 0.6777)
    omega_m = cosmo_params.get('omega_m', cosmo_params.get('omega_0', 0.3071))
    omega_wdm = cosmo_params.get('omega_wdm', omega_m)
    mu = 1.12
    alpha = 0.049 * (m_wdm_keV / 1.0)**(-1.11) * (omega_wdm / 0.25)**(0.11) * (h / 0.7)**(1.22)
    T_hm = 0.5
    k_hm = (1.0 / alpha) * (T_hm**(-mu / 5.0) - 1.0)**(1.0 / (2.0 * mu))
    rho_crit = 2.775e11
    rho_m = rho_crit * omega_m
    R_hm = np.pi / k_hm
    M_hm = (4.0 * np.pi / 3.0) * rho_m * (R_hm**3)
    return M_hm, k_hm

# =========================================================
# Main processing function
# =========================================================
def process_halo_data(input_file, cosmology_name, redshift, halo_def, relation):
    try:
        data = np.loadtxt(input_file, skiprows=1)
    except Exception as e:
        print(f"Error loading {input_file}: {e}")
        return None, None, None

    if data.ndim == 1:
        data = data.reshape(1, -1)

    Mvir = data[:, 0]
    z_data = np.full_like(data[:, 1], redshift, dtype=float)
    D_z = data[:, 2]
    sigma_M = data[:, 3]
    sigma_M_f = data[:, 4]
    sigma_M_f2 = data[:, 5]

    M_hm = 0.0
    if cosmology_name in cosmology_configs:
        cosmo_params = cosmology_configs[cosmology_name]
        if cosmo_params.get('cosmology_type') == 'LWDM':
            m_keV = cosmo_params.get('m_keV', 1.5)
            M_hm, _ = calculate_lwdm_half_mode_mass(cosmo_params, m_wdm_keV=m_keV)


    results = np.zeros((len(Mvir), 12))
    params = MODEL_PARAMS[halo_def][relation]

    for i in range(len(Mvir)):
        Mi = Mvir[i]
        zi = z_data[i]
        Dz = D_z[i]
        sig = sigma_M[i]
        sig_f = sigma_M_f[i]
        sig_f2 = sigma_M_f2[i]

        D_zf_sim = universal_D_zf(Dz, Mi, params, M_hm=M_hm)
        D_zf_eps = D_zf_eps_func(Dz, sig, sig_f)
        D_zf_eps2 = D_zf_eps_func(Dz, sig, sig_f2)

        c_sim = universal_concentration(Dz, D_zf_sim, sig, Mi, params, M_hm=M_hm)
        c_eps = universal_concentration(Dz, D_zf_eps, sig, Mi, params, M_hm=M_hm)
        c_eps2 = universal_concentration(Dz, D_zf_eps2, sig, Mi, params, M_hm=M_hm)

        results[i] = [
            Mi,          # 0
            sig,         # 1
            sig_f,       # 2
            zi,          # 3
            Dz,          # 4
            D_zf_sim,    # 5
            D_zf_eps,    # 6
            D_zf_eps2,   # 7
            c_sim,       # 8
            c_eps,       # 9
            c_eps2,      # 10
            M_hm         # 11
        ]

    return Mvir, D_z, results

# =========================================================
# Main
# =========================================================
def main():
    if len(sys.argv) < 6:
        sys.exit(1)

    input_file = sys.argv[1]
    cosmology = sys.argv[2]
    redshift = float(sys.argv[3])
    halo_def = sys.argv[4].lower()
    relation = sys.argv[5].lower()

    output_file = f"./output_halo_properties/halo_properties_{redshift}_{cosmology}_{halo_def}_{relation}.dat"

    Mvir, D_z, results = process_halo_data(input_file, cosmology, redshift, halo_def, relation)

    if results is None:
        sys.exit(1)

    header = (
        "# 0:Mvir[Msun/h] "
        "1:sigma(M,0) "
        "2:sigma(M_f,0)(f=0.50) "
        "3:z "
        "4:D(z) "
        "5:D(z_f)_simfit "
        "6:D(z_f)_eps(f=0.50) "
        "7:D(z_f)_eps(f=0.16) "
        "8:c_simfit "
        "9:c_eps(f=0.50) "
        "10:c_eps(f=0.16) "
        "11:M_hm[Msun/h]"
    )

    np.savetxt(
        output_file,
        results,
        fmt='%.6e %.6f %.6f %.3f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6e',
        header=header,
        comments=''
    )

if __name__ == '__main__':
    main()
