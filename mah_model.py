"""
Paper III: Universal Mass Accretion History (MAH) Model
=======================================================

Implements the parameter-free universal MAH model of Paper III and its coupling
with the individual-halo concentration relation of Paper II. Three tasks:

1. mah          : predict M(z)/M0
2. c_evolution  : predict the concentration history c(z)
3. reconstruct  : Bayesian MAP reconstruction of the MAH from (M_obs, c_obs, z_obs)

Formation-time / MAH source used by tasks 1 and 2:
  - merger tree given (--tree)      -> simulated data are used
  - input has 2 columns (Mass z_f)  -> simulated z_f is used
  - input has 1 column  (Mass)      -> omega_f is drawn from Paper II P(omega_f)
                                       and the universal MAH model is used

Outputs are written to ./output_mah_properties/ (one row per halo).
"""

import argparse
import contextlib
import io
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaincinv

from cosmology_params import cosmology_configs
from GrowthFactor import linear_growth
from VarianceCalculator import HaloMassVarianceCalculator
from halo_model import (
    MODEL_PARAMS,
    INDIVIDUAL_PARAMS,
    OMEGA_F_PDF_PARAMS,
    universal_D_zf,
    universal_concentration,
    calculate_lwdm_half_mode_mass,
    individual_concentration,
    sample_omega_f_sim,
    sample_omega_f_eps,
)

DELTA_SC0 = 1.686
LOG10_2 = np.log10(2.0)

MAH_PARAMS = {
    'vir':  {'omega_f_peak': {'CDM': 1.25, 'LWDM': 1.21}, 'a': 3.98},
    '200c': {'omega_f_peak': {'CDM': 1.29, 'LWDM': 1.27}, 'a': 3.06},
}

# Generalized Gamma PDF parameters for omega_f (Paper II Eq. 6) across vir and 200c
OMEGA_F_PDF_PARAMS = {
    'vir': {
        'CDM':  {'omega_0': 1.25, 'k': 2.14, 'p': 2.04},
        'LWDM': {'omega_0': 1.21, 'k': 1.82, 'p': 2.70}
    },
    '200c': {
        'CDM':  {'omega_0': 1.29, 'k': 2.14, 'p': 2.04},
        'LWDM': {'omega_0': 1.27, 'k': 1.82, 'p': 2.70}
    }
}

# Gaussian prior on u = log10(1 + omega_f) for Bayesian MAP inversion (Paper III Eq. 13)
# Under vir: mode=1.25, median=1.33 -> mu=log10(1+1.33)=0.368, sigma=0.094
# Under 200c: mode=1.29, median=1.37 -> mu=log10(1+1.37)=0.375, sigma=0.094
BAYES_PRIOR = {
    'vir':  {'mu': 0.368, 'sigma': 0.094},
    '200c': {'mu': 0.375, 'sigma': 0.094},
}

DEFAULT_Z_OUTPUTS = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
OUTPUT_DIR = './output_mah_properties'


# =========================================================
# Core Physics Functions (Paper III)
# =========================================================
def _pdf_key(cosmology_type):
    return 'LWDM' if cosmology_type == 'LWDM' else 'CDM'


def omega_f_peak(halo_def='vir', cosmology_type='CDM'):
    return MAH_PARAMS[halo_def]['omega_f_peak'][_pdf_key(cosmology_type)]


def omega_f_median(halo_def='vir', cosmology_type='CDM'):
    p = OMEGA_F_PDF_PARAMS[halo_def][_pdf_key(cosmology_type)]
    w0, k, pp = p['omega_0'], p['k'], p['p']
    return w0 * ((pp / k) * gammaincinv((k + 1.0) / pp, 0.5)) ** (1.0 / pp)


def A_hm_factor(M, M_hm=0.0, a=0.0):
    return 1.0 + a * np.asarray(M_hm, dtype=float) / np.asarray(M, dtype=float)


def beta_parameter(omega_f, M=1.0, M_hm=0.0, halo_def='vir', cosmology_type='CDM'):
    A_hm = A_hm_factor(M, M_hm, MAH_PARAMS[halo_def]['a'])
    return A_hm * np.asarray(omega_f, dtype=float) - omega_f_peak(halo_def, cosmology_type)


def tau_coordinate(D_z, D_zf, D_obs=1.0):
    return (1.0 / np.asarray(D_z) - 1.0 / D_obs) / (1.0 / np.asarray(D_zf) - 1.0 / D_obs)


def mah_ratio(tau, beta):
    tau = np.maximum(np.asarray(tau, dtype=float), 0.0)
    beta = np.asarray(beta, dtype=float)
    if beta.ndim > 0 and tau.ndim > beta.ndim:
        beta = beta[..., None]
    log_ratio = beta * np.log10(1.0 + tau) - (1.0 + beta) * tau * LOG10_2
    return 10.0 ** np.minimum(log_ratio, 0.0)


def D_zf_from_omega_f(omega_f, delta_S, D_obs=1.0):
    return 1.0 / (1.0 / D_obs + np.asarray(omega_f, dtype=float) * np.asarray(delta_S, dtype=float) / DELTA_SC0)


def omega_f_from_D(D_zf, delta_S, D_obs=1.0):
    return DELTA_SC0 * (1.0 / np.asarray(D_zf, dtype=float) - 1.0 / D_obs) / np.asarray(delta_S, dtype=float)


def bayesian_map_omega_f(y, alpha_omega, sigma_res, sigma_logc_obs=0.0,
                         mu_prior=BAYES_PRIOR['vir']['mu'], sigma_prior=BAYES_PRIOR['vir']['sigma']):
    var_like = sigma_res ** 2 + np.asarray(sigma_logc_obs, dtype=float) ** 2
    w_prior = 1.0 / sigma_prior ** 2
    w_like = alpha_omega ** 2 / var_like
    mu_post = (w_prior * mu_prior + w_like * np.asarray(y, dtype=float) / alpha_omega) / (w_prior + w_like)
    sigma_post = 1.0 / np.sqrt(w_prior + w_like)
    return mu_post, sigma_post


# =========================================================
# Unified Model Class (Papers I, II, & III)
# =========================================================
class UniversalHaloModel:
    def __init__(self, cosmology='planck18', halo_def='vir', verbose=False):
        if cosmology not in cosmology_configs:
            raise ValueError(f"Unknown cosmology '{cosmology}'. Available: {list(cosmology_configs.keys())}")
        if halo_def not in MAH_PARAMS:
            raise ValueError("halo_def must be 'vir' or '200c'")
        self.cosmology = cosmology
        self.halo_def = halo_def
        self.cosmo = cosmology_configs[cosmology]
        self.cosmology_type = self.cosmo['cosmology_type']
        self.verbose = verbose

        self.M_hm = 0.0
        if self.cosmology_type == 'LWDM':
            self.M_hm, _ = calculate_lwdm_half_mode_mass(self.cosmo, m_wdm_keV=self.cosmo.get('m_keV', 1.5))

        self._setup_growth()
        self._var = None
        self._M_lo = None
        self._M_hi = None

    def _setup_growth(self, z_max=200.0, n_grid=600):
        lz = np.linspace(0.0, np.log10(1.0 + z_max), n_grid)
        D, _ = linear_growth(
            10.0 ** lz - 1.0,
            h=self.cosmo['h'],
            omega_0=self.cosmo['omega_m'],
            w_0=self.cosmo['w_0'],
            w_a=self.cosmo['w_a'],
            cosmology=self.cosmology_type,
        )
        self._lz = lz
        self._lD = np.log10(D)

    def D(self, z):
        lz = np.log10(1.0 + np.asarray(z, dtype=float))
        lD = np.interp(lz, self._lz, self._lD)
        lD = np.where(lz > self._lz[-1], self._lD[-1] - (lz - self._lz[-1]), lD)
        return 10.0 ** lD

    def z_of_D(self, D):
        lD = np.log10(np.asarray(D, dtype=float))
        lz = np.interp(-lD, -self._lD, self._lz)
        lz = np.where(lD < self._lD[-1], self._lz[-1] + (self._lD[-1] - lD), lz)
        return 10.0 ** lz - 1.0

    def _ensure_sigma(self, M):
        M = np.asarray(M, dtype=float)
        M = M[np.isfinite(M) & (M > 0)]
        if M.size == 0:
            return
        lo, hi = M.min(), M.max()
        if self._var is not None and lo >= self._M_lo and hi <= self._M_hi:
            return
        M_lo = lo * 1e-5 if self._M_lo is None else min(self._M_lo, lo * 1e-5)
        M_hi = hi * 1e2 if self._M_hi is None else max(self._M_hi, hi * 1e2)
        sink = sys.stdout if self.verbose else io.StringIO()
        with contextlib.redirect_stdout(sink):
            var = HaloMassVarianceCalculator(
                cosmology_name=self.cosmology, auto_adjust=True, safety_factor=0.1, table_size=1000
            )
            var.initialize_for_masses(np.array([M_lo, M_hi]))
        self._var, self._M_lo, self._M_hi = var, M_lo, M_hi

    def sigma(self, M):
        self._ensure_sigma(M)
        R = self._var._calculate_radius_from_mass(np.asarray(M, dtype=float))
        return np.sqrt(self._var.sigma2(R))

    def delta_S(self, M):
        M = np.asarray(M, dtype=float)
        return np.sqrt(np.maximum(self.sigma(0.5 * M) ** 2 - self.sigma(M) ** 2, 1e-30))

    # --- Paper I & II ---
    def concentration_paper1(self, M, z=0.0, statistic='peak'):
        params = MODEL_PARAMS[self.halo_def][statistic]
        M = np.asarray(M, dtype=float)
        D_z = self.D(z)
        D_zf = universal_D_zf(D_z, M, params, M_hm=self.M_hm)
        return universal_concentration(D_z, D_zf, self.sigma(M), M, params, M_hm=self.M_hm)

    def concentration_paper2(self, M, omega_f, z=0.0, relation='relaxed'):
        params = INDIVIDUAL_PARAMS[self.halo_def][relation]
        M = np.asarray(M, dtype=float)
        return individual_concentration(self.D(z), self.sigma(M), omega_f, M, params, M_hm=self.M_hm)

    def omega_f(self, M, z_f, z=0.0):
        return omega_f_from_D(self.D(z_f), self.delta_S(M), self.D(z))

    def formation_redshift(self, M, omega_f, z=0.0):
        D_zf = D_zf_from_omega_f(np.maximum(omega_f, 1e-3), self.delta_S(M), self.D(z))
        return self.z_of_D(D_zf)

    def sample_omega_f(self, n, method='sim', seed=42):
        rng = np.random.default_rng(seed)
        if method in ['sim', 'mc_sim']:
            p_pdf = OMEGA_F_PDF_PARAMS[self.halo_def][_pdf_key(self.cosmology_type)]
            w0, k, p = p_pdf['omega_0'], p_pdf['k'], p_pdf['p']
            u = rng.gamma(shape=(k + 1.0) / p, scale=1.0, size=n)
            return w0 * ((p / k) * u) ** (1.0 / p)
        if method in ['eps', 'mc_eps']:
            return sample_omega_f_eps(n, rng=rng)
        if method == 'median':
            return np.full(n, omega_f_median(self.halo_def, self.cosmology_type))
        raise ValueError(f"Unknown omega_f sampling method '{method}'")

    # --- Paper III: Universal MAH ---
    def beta(self, M, omega_f):
        return beta_parameter(omega_f, M, self.M_hm, self.halo_def, self.cosmology_type)

    def mah(self, z, M, omega_f=None, z_f=None, z_obs=0.0):
        """Evaluate M(z)/M(z_obs) for given halo masses across redshift array z."""
        z = np.atleast_1d(np.asarray(z, dtype=float))
        M = np.atleast_1d(np.asarray(M, dtype=float))
        D_obs = float(self.D(z_obs))
        dS = self.delta_S(M)
        if z_f is not None:
            z_f = np.broadcast_to(np.asarray(z_f, dtype=float), M.shape)
            omega = omega_f_from_D(self.D(z_f), dS, D_obs)
        else:
            if omega_f is None:
                omega_f = omega_f_median(self.halo_def, self.cosmology_type)
            omega = np.broadcast_to(np.asarray(omega_f, dtype=float), M.shape).copy()
        omega_safe = np.maximum(omega, 1e-3)
        D_zf = D_zf_from_omega_f(omega_safe, dS, D_obs)
        beta = self.beta(M, omega_safe)

        tau = tau_coordinate(self.D(z)[None, :], D_zf[:, None], D_obs)
        ratio = mah_ratio(tau, beta)
        ratio[:, z < z_obs - 1e-10] = np.nan
        return ratio

    # --- Paper III Sec 4: Forward c(z) ---
    def concentration_track(self, z_eval, z_mah, M_mah, relation='relaxed'):
        params = INDIVIDUAL_PARAMS[self.halo_def][relation]
        z_eval = np.atleast_1d(np.asarray(z_eval, dtype=float))
        z_m = np.asarray(z_mah, dtype=float)
        M_m = np.asarray(M_mah, dtype=float)
        ok = np.isfinite(z_m) & np.isfinite(M_m) & (M_m > 0)
        order = np.argsort(z_m[ok])
        z_m, M_m = z_m[ok][order], M_m[ok][order]
        nan = np.full(z_eval.shape, np.nan)
        if z_m.size < 2:
            return nan

        lM = np.log10(np.minimum.accumulate(M_m))
        lM_e = np.interp(z_eval, z_m, lM, left=np.nan, right=np.nan)
        z_fp = np.interp(lM_e - LOG10_2, lM[::-1], z_m[::-1], left=np.nan, right=np.nan)

        good = np.isfinite(lM_e) & np.isfinite(z_fp)
        c = nan.copy()
        if np.any(good):
            Mg = 10.0 ** lM_e[good]
            omega = np.maximum(omega_f_from_D(self.D(z_fp[good]), self.delta_S(Mg), self.D(z_eval[good])), 0.0)
            c[good] = individual_concentration(self.D(z_eval[good]), self.sigma(Mg), omega, Mg, params, M_hm=self.M_hm)
        return c

    def concentration_track_model(self, z_eval, M, omega_f=None, z_f=None, z_obs=0.0, relation='relaxed'):
        """Predict continuous concentration c(z) along the progenitor branch using the universal MAH model (vectorized)."""
        M = np.atleast_1d(np.asarray(M, dtype=float))
        z_eval = np.atleast_1d(np.asarray(z_eval, dtype=float))
        D_obs = float(self.D(z_obs))
        dS = self.delta_S(M)

        if z_f is not None:
            z_f = np.broadcast_to(np.asarray(z_f, dtype=float), M.shape)
            omega = omega_f_from_D(self.D(z_f), dS, D_obs)
        else:
            if omega_f is None:
                omega_f = omega_f_median(self.halo_def, self.cosmology_type)
            omega = np.broadcast_to(np.asarray(omega_f, dtype=float), M.shape).copy()

        omega_safe = np.maximum(omega, 1e-3)
        D_zf = D_zf_from_omega_f(omega_safe, dS, D_obs)
        beta = self.beta(M, omega_safe)

        D_z = self.D(z_eval)
        tau = tau_coordinate(D_z[None, :], D_zf[:, None], D_obs)
        log_ratio = beta[:, None] * np.log10(1.0 + tau) - (1.0 + beta[:, None]) * tau * LOG10_2
        M_z = M[:, None] * (10.0 ** np.minimum(log_ratio, 0.0))

        y_target = log_ratio - LOG10_2
        ln10 = np.log(10.0)
        tau_h = tau + 1.0
        b_2d = beta[:, None]
        for _ in range(6):
            tau_h = np.maximum(tau_h, 0.0)
            f_val = b_2d * np.log10(1.0 + tau_h) - (1.0 + b_2d) * tau_h * LOG10_2
            df_val = b_2d / ((1.0 + tau_h) * ln10) - (1.0 + b_2d) * LOG10_2
            tau_h = tau_h - (f_val - y_target) / df_val
        tau_h = np.maximum(tau_h, 0.0)

        inv_D_zf_prog = 1.0 / D_obs + tau_h * (1.0 / D_zf[:, None] - 1.0 / D_obs)
        D_zf_prog = 1.0 / inv_D_zf_prog

        dS_z = self.delta_S(M_z)
        omega_z = np.maximum(DELTA_SC0 * (1.0 / D_zf_prog - 1.0 / D_z[None, :]) / dS_z, 0.0)

        params = INDIVIDUAL_PARAMS[self.halo_def][relation]
        sigma_z = self.sigma(M_z)
        D_z_broad = np.broadcast_to(D_z[None, :], M_z.shape)
        c = individual_concentration(D_z_broad, sigma_z, omega_z, M_z, params, M_hm=self.M_hm)
        c[:, z_eval < z_obs - 1e-10] = np.nan
        return c

    # --- Paper III Sec 5: Bayesian Inversion ---
    def _solve_M0(self, M_obs, omega_f, z_obs):
        if z_obs < 1e-6:
            return float(M_obs)
        l_obs = np.log10(M_obs)
        def f(lM0):
            M0 = 10.0 ** lM0
            D_zf = D_zf_from_omega_f(omega_f, self.delta_S(M0), 1.0)
            tau = tau_coordinate(self.D(z_obs), D_zf, 1.0)
            ratio = mah_ratio(tau, self.beta(M0, omega_f))
            return lM0 + np.log10(float(ratio)) - l_obs

        a, b = l_obs, l_obs + 0.5
        if f(a) >= 0.0:
            return float(M_obs)
        while f(b) <= 0.0 and b < l_obs + 4.0:
            b += 0.5
        return float(10.0 ** brentq(f, a, b, xtol=1e-6))

    def reconstruct(self, M_obs, c_obs, z_obs=0.0, z_eval=None, relation='relaxed', sigma_logc_obs=0.0):
        """Reconstruct the most probable (MAP) MAH and z_f from (M_obs, c_obs, z_obs).
        sigma_logc_obs: optional observational error of log10(c_obs), added in quadrature."""
        M_obs = np.atleast_1d(np.asarray(M_obs, dtype=float))
        c_obs = np.broadcast_to(np.asarray(c_obs, dtype=float), M_obs.shape)
        z_obs = np.broadcast_to(np.asarray(z_obs, dtype=float), M_obs.shape)
        p = INDIVIDUAL_PARAMS[self.halo_def][relation]

        y = (np.log10(c_obs) - p['A'] - p['alpha_D'] * np.log10(self.D(z_obs))
             - p['alpha_sigma'] * np.log10(self.sigma(M_obs))
             - p['gamma'] * np.log10(A_hm_factor(M_obs, self.M_hm, p['a'])))
        prior = BAYES_PRIOR[self.halo_def]
        mu_post, _ = bayesian_map_omega_f(y, p['beta'], p['sigma_res'], sigma_logc_obs=sigma_logc_obs,
                                          mu_prior=prior['mu'], sigma_prior=prior['sigma'])
        omega_map = np.maximum(10.0 ** mu_post - 1.0, 0.01)

        M0_map = np.array([self._solve_M0(M_obs[i], omega_map[i], z_obs[i]) for i in range(len(M_obs))])
        zf_map = self.formation_redshift(M0_map, omega_map, z=0.0)

        out = {
            'omega_f_MAP': omega_map,
            'z_f_MAP': zf_map,
            'M0_MAP': M0_map
        }
        if z_eval is not None:
            z_eval = np.atleast_1d(np.asarray(z_eval, dtype=float))
            tracks = np.full((len(M_obs), len(z_eval)), np.nan)
            for i in range(len(M_obs)):
                D_zf = D_zf_from_omega_f(omega_map[i], self.delta_S(M0_map[i]), 1.0)
                tau = tau_coordinate(self.D(z_eval), D_zf, 1.0)
                tracks[i] = (M0_map[i] / M_obs[i]) * mah_ratio(tau, self.beta(M0_map[i], omega_map[i]))
            out['tracks'] = tracks
        return out



# =========================================================
# File Loaders
# =========================================================
def _load_input(path):
    if not os.path.exists(path):
        sys.exit(f"ERROR: Input file '{path}' not found.")
    return np.loadtxt(path, ndmin=2, comments='#')


def _parse_z_list(z_str):
    if z_str is None:
        return list(DEFAULT_Z_OUTPUTS)
    return [float(x.strip()) for x in z_str.split(',') if x.strip()]


def _load_tree(path):
    """Read a merger-tree file (columns: halo_id  z  M(z)) and return a list of
    (halo_id, z_array, M_array) sorted by increasing redshift."""
    tree = _load_input(path)
    if tree.shape[1] < 3:
        sys.exit("ERROR: Merger-tree file must have 3 columns: halo_id  z  M(z)")
    halos = []
    for hid in np.unique(tree[:, 0]):
        sel = tree[:, 0] == hid
        order = np.argsort(tree[sel, 1])
        halos.append((int(hid), tree[sel, 1][order], tree[sel, 2][order]))
    return halos


def _zf_from_track(z, M):
    """Half-mass formation redshift measured directly from a simulated main branch."""
    lM = np.log10(np.minimum.accumulate(M))           # enforce monotonic main branch
    target = lM[0] - LOG10_2
    if lM[-1] > target:
        return np.nan                                  # tree too short to reach M0/2
    return float(np.interp(target, lM[::-1], z[::-1]))


# =========================================================
# Task 1: MAH prediction  M(z)/M0
# =========================================================
def run_task_mah(args):
    """Formation-time source (in order of priority):
       1. --tree file           -> z_f measured from the simulated merger tree
       2. input with 2 columns  -> z_f given by the user (simulation)
       3. input with 1 column   -> omega_f drawn from Paper II P(omega_f)   [default]
          (or --mode eps / median)"""
    model = UniversalHaloModel(args.cosmology, args.halo_def)
    z_obs = float(args.z)
    z_eval = _parse_z_list(args.z_list)

    if args.tree:
        halos = _load_tree(args.tree)
        M = np.array([m[0] for _, _, m in halos])
        z_f = np.array([_zf_from_track(z, m) for _, z, m in halos])
        omega_f = model.omega_f(M, z_f, z=z_obs)
        label = "tree"
    else:
        data = _load_input(args.input)
        M = data[:, 0]
        mode = args.mode
        if mode == 'auto':
            mode = 'sim_zf' if data.shape[1] >= 2 else 'sim'
        if mode == 'sim_zf':
            if data.shape[1] < 2:
                sys.exit("ERROR: --mode sim_zf requires a 2-column input (Mass  z_f).")
            z_f = data[:, 1]
            omega_f = model.omega_f(M, z_f, z=z_obs)
            label = "sim_zf"
        elif mode in ('sim', 'eps', 'median'):
            omega_f = model.sample_omega_f(len(M), method=mode, seed=args.seed)
            z_f = model.formation_redshift(M, omega_f, z=z_obs)
            label = {'sim': 'mc_sim', 'eps': 'mc_eps', 'median': 'median'}[mode]
        else:
            sys.exit(f"ERROR: Unrecognized mode '{args.mode}'")

    beta = model.beta(M, omega_f)
    ratios = model.mah(z_eval, M, omega_f=omega_f, z_obs=z_obs)

    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir, f"mah_{z_obs}_{args.cosmology}_{args.halo_def}_{label}.dat")
    header = " ".join(["0:M0[Msun/h]", "1:z_f", "2:omega_f", "3:beta"]
                      + [f"{4 + i}:M(z={zv:g})/M0" for i, zv in enumerate(z_eval)])
    table = np.column_stack([M, z_f, omega_f, beta, ratios])
    fmt = ['%.6e', '%.4f', '%.4f', '%.4f'] + ['%.5e'] * len(z_eval)
    np.savetxt(out_file, table, fmt=fmt, header=header, comments='# ')
    print(f"[Paper III] MAH ({label}) saved to: {out_file}")
    return out_file


# =========================================================
# Task 2: Forward concentration evolution c(z)
# =========================================================
def run_task_c_evolution(args):
    """MAH source (in order of priority):
       1. --tree file           -> c(z) computed along the simulated main branch
       2. input with 2 columns  -> universal MAH model with the user-given z_f
       3. input with 1 column   -> universal MAH model with omega_f from Paper II P(omega_f) [default]"""
    model = UniversalHaloModel(args.cosmology, args.halo_def)
    z_obs = float(args.z)
    z_eval = _parse_z_list(args.z_list)
    if z_obs not in z_eval:
        z_eval = sorted([z_obs] + z_eval)

    if args.tree:
        halos = _load_tree(args.tree)
        M = np.array([m[0] for _, _, m in halos])
        z_f = np.array([_zf_from_track(z, m) for _, z, m in halos])
        c_out = np.array([model.concentration_track(z_eval, z, m, relation=args.relation)
                          for _, z, m in halos])
        label = "tree"
    else:
        data = _load_input(args.input)
        M = data[:, 0]
        if data.shape[1] >= 2:
            z_f = data[:, 1]
            omega_f = model.omega_f(M, z_f, z=z_obs)
            label = "sim_zf"
        else:
            omega_f = model.sample_omega_f(len(M), method='sim', seed=args.seed)
            z_f = model.formation_redshift(M, omega_f, z=z_obs)
            label = "mc_sim"
        c_out = model.concentration_track_model(z_eval, M, omega_f=omega_f, z_obs=z_obs,
                                                relation=args.relation)

    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir,
                            f"c_evolution_{z_obs}_{args.cosmology}_{args.halo_def}_{args.relation}_{label}.dat")
    header = " ".join(["0:M0[Msun/h]", "1:z_f"]
                      + [f"{2 + i}:c(z={zv:g})" for i, zv in enumerate(z_eval)])
    table = np.column_stack([M, z_f, c_out])
    fmt = ['%.6e', '%.4f'] + ['%.3f'] * len(z_eval)
    np.savetxt(out_file, table, fmt=fmt, header=header, comments='# ')
    print(f"[Paper III] c(z) evolution ({label}) saved to: {out_file}")
    return out_file


# =========================================================
# Task 3: Bayesian MAP reconstruction of the MAH
# =========================================================
def run_task_reconstruct(args):
    """Input columns: M_obs  c_obs  [z_obs]  [sigma_log10_c_obs]"""
    model = UniversalHaloModel(args.cosmology, args.halo_def)
    data = _load_input(args.input)
    if data.shape[1] < 2:
        sys.exit("ERROR: Reconstruction input needs at least 2 columns: M_obs  c_obs  [z_obs]  [sigma_logc]")
    M_obs, c_obs = data[:, 0], data[:, 1]
    z_obs = data[:, 2] if data.shape[1] >= 3 else np.full_like(M_obs, float(args.z))
    sig_c = data[:, 3] if data.shape[1] >= 4 else 0.0
    z_eval = _parse_z_list(args.z_list)

    recon = model.reconstruct(M_obs, c_obs, z_obs=z_obs, z_eval=z_eval,
                              relation=args.relation, sigma_logc_obs=sig_c)

    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir,
                            f"reconstruct_{args.cosmology}_{args.halo_def}_{args.relation}.dat")
    header = " ".join(["0:M_obs[Msun/h]", "1:c_obs", "2:z_obs", "3:z_f_MAP", "4:omega_f_MAP"]
                      + [f"{5 + i}:M(z={zv:g})/M_obs" for i, zv in enumerate(z_eval)])
    table = np.column_stack([M_obs, c_obs, z_obs, recon['z_f_MAP'], recon['omega_f_MAP'], recon['tracks']])
    fmt = ['%.6e', '%.3f', '%.3f', '%.4f', '%.4f'] + ['%.5e'] * len(z_eval)
    np.savetxt(out_file, table, fmt=fmt, header=header, comments='# ')
    print(f"[Paper III] Reconstructed MAHs saved to: {out_file}")
    return out_file


# =========================================================
# Command-line interface
# =========================================================
def main():
    parser = argparse.ArgumentParser(description="Paper III: universal MAH, forward c(z), and Bayesian MAH reconstruction.")
    sub = parser.add_subparsers(dest='task', required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--cosmology', default='planck18', help="Cosmology in cosmology_params.py (default: planck18)")
    common.add_argument('--halo-def', dest='halo_def', default='vir', choices=['vir', '200c'], help="Halo definition (default: vir)")
    common.add_argument('--z', type=float, default=0.0, help="Redshift of the input halos (default: 0)")
    common.add_argument('--z-list', dest='z_list', default=None, help="Output redshifts, e.g. '0.5,1,2,3,4,5'")
    common.add_argument('--outdir', default=OUTPUT_DIR, help=f"Output directory (default: {OUTPUT_DIR})")
    common.add_argument('--seed', type=int, default=42, help="Random seed for Monte Carlo sampling (default: 42)")

    rel_kw = dict(default='relaxed', choices=['relaxed', 'all_halos', 'relaxed_einasto'],
                  help="Paper II concentration relation (default: relaxed)")

    p1 = sub.add_parser('mah', parents=[common], help="Predict M(z)/M0")
    p1.add_argument('--input', default='input_masses.dat', help="Columns: Mass [z_f]")
    p1.add_argument('--tree', default=None, help="Optional merger tree (columns: halo_id z M(z)); z_f is measured from it")
    p1.add_argument('--mode', default='auto', choices=['auto', 'sim_zf', 'sim', 'eps', 'median'],
                    help="auto: use z_f if given, else sample Paper II P(omega_f); eps: EPS P(omega_f); median: median omega_f")
    p1.set_defaults(func=run_task_mah)

    p2 = sub.add_parser('c_evolution', parents=[common], help="Predict c(z) along the main branch")
    p2.add_argument('--input', default='input_masses.dat', help="Columns: Mass [z_f]")
    p2.add_argument('--tree', default=None, help="Optional merger tree (columns: halo_id z M(z))")
    p2.add_argument('--relation', **rel_kw)
    p2.set_defaults(func=run_task_c_evolution)

    p3 = sub.add_parser('reconstruct', parents=[common], help="Reconstruct the MAH from (M_obs, c_obs)")
    p3.add_argument('--input', default='input_mass_concentration.dat', help="Columns: M_obs c_obs [z_obs] [sigma_logc]")
    p3.add_argument('--relation', **rel_kw)
    p3.set_defaults(func=run_task_reconstruct)

    args = parser.parse_args()
    if getattr(args, 'relation', None) == 'relaxed_einasto' and args.halo_def != '200c':
        sys.exit("ERROR: 'relaxed_einasto' is calibrated for 200c only (Paper II Eq. 8).")
    args.func(args)


if __name__ == '__main__':
    main()
