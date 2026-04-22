# """
# Halo Mass Variance Calculator
# ==============================

# This module provides a class for computing the variance of halo mass fluctuations
# based on cosmological power spectra. It supports multiple cosmological models
# including LCDM, SCDM, OCDM, LWDM, wCDM, and Scale-Free cosmologies.

# Features:
# ---------
# - Support for multiple cosmological models via configuration system
# - Integration with CAMB for accurate power spectrum computation
# - Interpolation tables for fast variance calculations
# - Automatic adjustment of interpolation ranges based on input masses
# - Vectorized operations for efficient batch processing
# """
 
import numpy as np
from scipy import integrate
import camb
from scipy.interpolate import interp1d
from GrowthFactor import linear_growth
from cosmology_params import cosmology_configs

class HaloMassVarianceCalculator:
    def __init__(self, cosmology_name='LCDM', table_size=1000, R_min=None, R_max=None,
                 auto_adjust=True, safety_factor=0.1,
                 window_type=None, sharpk_c=2.42,
                 tophat_nk=2048, tophat_tail_factor=200.0):
        """Initialize with a cosmology configuration name and interpolation table"""
        if cosmology_name not in cosmology_configs:
            raise ValueError(
                f"Unknown cosmology name: {cosmology_name}. "
                f"Available options are: {list(cosmology_configs.keys())}"
            )

        self.cosmo = cosmology_configs[cosmology_name]
        self.cosmo_name = cosmology_name

        # Lookup table parameters
        self.table_size = table_size
        self.auto_adjust = auto_adjust
        self.safety_factor = safety_factor

        # Window-function control
        self.sharpk_c = sharpk_c
        if window_type is None:
            self.window_type = 'sharpk' if self.cosmo['cosmology_type'] == 'LWDM' else 'tophat'
        else:
            if window_type not in ['tophat', 'sharpk']:
                raise ValueError("window_type must be None, 'tophat', or 'sharpk'")
            self.window_type = window_type

        # Fast top-hat integration control for CAMB-based cosmologies
        self.tophat_nk = int(tophat_nk)
        self.tophat_tail_factor = float(tophat_tail_factor)

        # Initialize table bounds
        self.R_min = R_min if R_min is not None else 0.001
        self.R_max = R_max if R_max is not None else 50.0
        self.R_table = np.logspace(np.log10(self.R_min), np.log10(self.R_max), table_size)

        self.sigma2_table = None
        self.sigma2_interp = None
        self.sigma_table_initialized = False
        self.AA_norm = None

        self.rho_c = 3.0e4 / (8.0 * np.pi * 6.67430e-11) * 1.5513826e-2
        self.rho_0 = self.rho_c * self.cosmo['omega_m']

        # Power spectrum related containers
        self.power_spectrum_func = None
        self.kmin_pk = None
        self.kmax_pk = None
        self.k_grid_tophat = None
        self.pk_grid_tophat = None

        # Compute power spectrum using CAMB if needed
        self._setup_power_spectrum()

    def _calculate_radius_from_mass(self, M):
        """Calculate filter radius from halo mass"""
        M = np.asarray(M, dtype=np.float64)

        # top-hat: M = (4/3) pi rho_0 R^3
        if self.window_type == 'tophat':
            return (3.0 * M / (4.0 * np.pi * self.rho_0)) ** (1.0 / 3.0)

        # sharp-k: M = (4/3) pi rho_0 (c R)^3
        elif self.window_type == 'sharpk':
            return (3.0 * M / (4.0 * np.pi * self.rho_0)) ** (1.0 / 3.0) / self.sharpk_c

        else:
            raise ValueError(f"Unknown window_type: {self.window_type}")

    def _adjust_table_bounds(self, masses):
        """Adjust R_min and R_max based on provided halo masses"""
        if not self.auto_adjust:
            return

        radii = self._calculate_radius_from_mass(np.array(masses))
        radii = radii[np.isfinite(radii) & (radii > 0)]

        if len(radii) == 0:
            return

        R_data_min = np.min(radii)
        R_data_max = np.max(radii)

        # Apply safety factor to extend bounds
        R_range = max(R_data_max - R_data_min, 1e-12)
        self.R_min = max(1e-6, R_data_min - self.safety_factor * R_range)
        self.R_max = R_data_max + self.safety_factor * R_range

        print(f"Auto-adjusted radius bounds: R_min={self.R_min:.4e}, R_max={self.R_max:.4e}")
        print(f"Data radius range: {R_data_min:.4e} - {R_data_max:.4e}")

        self.R_table = np.logspace(np.log10(self.R_min), np.log10(self.R_max), self.table_size)

    def _setup_power_spectrum(self):
        """Setup power spectrum using CAMB"""
        params = self.cosmo

        # For scale-free cosmologies
        if params['cosmology_type'] == 'EdS':
            self.power_spectrum_func = None
            return

        # For other cosmologies, use CAMB
        pars = camb.CAMBparams()

        if params['cosmology_type'] in ['CDM', 'OCDM', 'wCDM', 'LWDM']:
            pars.set_cosmology(
                H0=params['h'] * 100,
                ombh2=params['omega_b'] * params['h']**2,
                omch2=(params['omega_m'] - params['omega_b']) * params['h']**2,
                omk=0 if params['cosmology_type'] in ['CDM', 'LWDM', 'wCDM'] else 1 - params['omega_m']
            )

            if params['cosmology_type'] == 'wCDM':
                pars.set_dark_energy(w=params['w_0'], wa=params['w_a'])

            pars.InitPower.set_params(ns=params['n'])
        else:
            raise ValueError(f"Unsupported cosmology type: {params['cosmology_type']}")

        # top-hat branch: extend k range a bit
        if self.window_type == 'tophat':
            camb_kmax = 200.0
            npoints = 2000
            minkh = 1e-5
        else:
            camb_kmax = 100.0
            npoints = 1000
            minkh = 1e-4

        pars.set_matter_power(
            redshifts=[0.0],
            kmax=camb_kmax,
            k_per_logint=0
        )

        results = camb.get_results(pars)
        kh, z, P_k = results.get_matter_power_spectrum(
            minkh=minkh,
            maxkh=camb_kmax,
            npoints=npoints
        )
        P_k = P_k[0]  # Take z=0

        # For LWDM model, apply WDM correction to power spectrum
        if params['cosmology_type'] == 'LWDM':
            m_wdm_keV = params.get('m_keV', 0.8)
            Omega_wdm = params.get('omega_wdm', params['omega_m'])
            h = params['h']

            # Viel+2005 formula
            nu = 1.12
            alpha = 0.049 * (m_wdm_keV / 1.0)**(-1.11) * (Omega_wdm / 0.25)**(0.11) * (h / 0.7)**(1.22)

            # WDM transfer function correction factor
            Tk = (1.0 + (alpha * kh)**(2.0 * nu))**(-5.0 / nu)
            P_k = P_k * Tk**2

        # Keep only valid data
        valid = np.isfinite(kh) & np.isfinite(P_k) & (kh > 0) & (P_k > 0)
        kh = kh[valid]
        P_k = P_k[valid]

        if len(kh) < 10:
            raise RuntimeError("Power spectrum grid is too small or invalid.")

        self.kmin_pk = kh.min()
        self.kmax_pk = kh.max()

        self.power_spectrum_func = interp1d(
            kh, P_k,
            kind='linear',
            bounds_error=False,
            fill_value=(P_k[0], P_k[-1])
        )

        # Precompute top-hat integration grid for fast CAMB-based calculations
        if self.window_type == 'tophat':
            self.k_grid_tophat = np.logspace(
                np.log10(self.kmin_pk),
                np.log10(self.kmax_pk),
                self.tophat_nk
            )
            self.pk_grid_tophat = self.power_spectrum_func(self.k_grid_tophat)

    def _window_function_squared(self, k, R):
        """Return |W(kR)|^2 for the current window type."""
        k = np.asarray(k, dtype=np.float64)
        kR = k * R

        if self.window_type == 'tophat':
            W = np.ones_like(kR, dtype=np.float64)
            mask = np.abs(kR) >= 1e-8
            kr = kR[mask]
            W[mask] = 3.0 * (np.sin(kr) - kr * np.cos(kr)) / (kr ** 3)
            W2 = W * W
            return float(W2) if W2.ndim == 0 else W2

        elif self.window_type == 'sharpk':
            W2 = np.where(k < 1.0 / R, 1.0, 0.0)
            return float(W2) if W2.ndim == 0 else W2

        else:
            raise ValueError(f"Unknown window_type: {self.window_type}")

    def _sigma2_integrand(self, k, R, AA):
        """Integrand function for sigma2 calculation"""
        W2 = self._window_function_squared(k, R)

        if np.isscalar(W2) and W2 == 0.0:
            return 0.0

        if self.cosmo['cosmology_type'] == 'EdS':
            return (1.0 / (2.0 * np.pi * np.pi)) * AA * k**(2 + self.cosmo['n']) * W2
        else:
            Pk_val = self.power_spectrum_func(k)
            return (1.0 / (2.0 * np.pi * np.pi)) * AA * k**2 * Pk_val * W2

    def _integrate_sigma2_tophat_fast(self, R, AA):
        """
        Fast and stable sigma^2(R) integration for top-hat window.
        Used only for CAMB-based cosmologies.
        """
        if self.cosmo['cosmology_type'] == 'EdS':
            raise RuntimeError("EdS cosmology should not use _integrate_sigma2_tophat_fast().")

        if self.kmin_pk is None or self.kmax_pk is None:
            raise RuntimeError("Power spectrum grid is not initialized for fast top-hat integration.")

        R = float(R)
        if R <= 0:
            return 0.0

        k_upper = min(self.kmax_pk, max(self.tophat_tail_factor / R, self.kmin_pk * 1.01))

        if k_upper <= self.kmin_pk:
            k = np.logspace(np.log10(self.kmin_pk), np.log10(self.kmin_pk * 1.01), 32)
            pk = self.power_spectrum_func(k)
        else:
            mask = self.k_grid_tophat <= k_upper
            k = self.k_grid_tophat[mask]
            pk = self.pk_grid_tophat[mask]

            if k.size < 64:
                k = np.logspace(np.log10(self.kmin_pk), np.log10(k_upper), 256)
                pk = self.power_spectrum_func(k)

        W2 = self._window_function_squared(k, R)
        integrand = (1.0 / (2.0 * np.pi * np.pi)) * AA * k**2 * pk * W2

        return integrate.simpson(integrand, x=k)

    def _integrate_sigma2(self, R, AA):
        """Numerically integrate sigma^2(R)"""
        R = float(R)
        if R <= 0:
            return 0.0

        # Scale-free EdS branch: always use direct numerical integration
        if self.cosmo['cosmology_type'] == 'EdS':
            if self.window_type == 'sharpk':
                k_upper = max(1.0 / R, 1e-8)
            else:
                # finite upper limit for numerical top-hat integration
                k_upper = max(self.tophat_tail_factor / R, 1e-8)

            result, _ = integrate.quad(
                self._sigma2_integrand,
                0.0, k_upper,
                args=(R, AA),
                epsabs=1.49e-6,
                epsrel=1.49e-6,
                limit=1000
            )
            return result

        # LWDM / other sharp-k branch
        if self.window_type == 'sharpk':
            k_upper = max(1.0 / R, 1e-8)
            result, _ = integrate.quad(
                self._sigma2_integrand,
                0.0, k_upper,
                args=(R, AA),
                epsabs=1.49e-6,
                epsrel=1.49e-6,
                limit=1000
            )
            return result

        # CAMB-based top-hat branch
        return self._integrate_sigma2_tophat_fast(R, AA)

    def _setup_sigma2_table(self):
        """Precompute sigma2 values for lookup table"""
        print("Precomputing sigma2 lookup table...")
        print(f"Table range: R_min={self.R_min:.4e}, R_max={self.R_max:.4e}, size={self.table_size}")
        print(f"Window type: {self.window_type}")

        # First, normalize to sigma_8
        AA_temp = 1.0
        s2_8 = self._integrate_sigma2(8.0, AA_temp)

        if not np.isfinite(s2_8) or s2_8 <= 0:
            raise RuntimeError(f"Invalid sigma^2(8): {s2_8}")

        self.AA_norm = self.cosmo['sigma_8'] ** 2 / s2_8

        # Precompute sigma2 for all R values in the table
        self.sigma2_table = np.zeros(self.table_size)
        for i, R in enumerate(self.R_table):
            result = self._integrate_sigma2(R, AA_temp)
            if not np.isfinite(result) or result < 0:
                raise RuntimeError(f"Invalid sigma^2 result at R={R:.4e}: {result}")
            self.sigma2_table[i] = result * self.AA_norm

        self.sigma2_interp = interp1d(
            self.R_table,
            self.sigma2_table,
            kind='cubic',
            bounds_error=False,
            fill_value=(self.sigma2_table[0], self.sigma2_table[-1])
        )

        self.sigma_table_initialized = True
        print("Sigma2 lookup table initialized successfully.")

    def initialize_for_masses(self, masses):
        """Initialize or reinitialize the lookup table for specific halo masses"""
        print("Initializing lookup table for provided halo masses...")

        if self.auto_adjust:
            self._adjust_table_bounds(masses)

        self._setup_sigma2_table()

    def sigma2(self, R, AA=None):
        """Calculate sigma^2(R) using interpolation lookup table"""
        if not self.sigma_table_initialized:
            raise RuntimeError("Sigma2 lookup table not initialized. Call initialize_for_masses first.")

        R = np.asarray(R, dtype=np.float64)
        R_clipped = np.clip(R, self.R_min, self.R_max)
        return self.sigma2_interp(R_clipped)

    def calculate_variance(self, M, redshift, return_all=True):
        M = np.asarray(M, dtype=np.float64)

        if not self.sigma_table_initialized:
            self.initialize_for_masses(M)

        Dz, _ = linear_growth(
            redshift,
            h=self.cosmo['h'],
            omega_0=self.cosmo['omega_m'],
            w_0=self.cosmo['w_0'],
            w_a=self.cosmo['w_a'],
            cosmology=self.cosmo['cosmology_type']
        )

        valid_mask = np.isfinite(M) & (M > 0)

        results = {
            'M': M,
            'D_z': np.full_like(M, np.nan),
            'sigma_M': np.full_like(M, np.nan),
            'sigma_M_f': np.full_like(M, np.nan),
            'sigma_M_f2': np.full_like(M, np.nan),
        }

        if np.any(valid_mask):
            valid_M = M[valid_mask]

            # Use the same mass-radius mapping as the chosen filter
            R = self._calculate_radius_from_mass(valid_M)
            Rf = self._calculate_radius_from_mass(0.5 * valid_M)
            Rf2 = self._calculate_radius_from_mass(0.16 * valid_M)

            sigma2_M = self.sigma2(R)
            sigma2_M_f = self.sigma2(Rf)
            sigma2_M_f2 = self.sigma2(Rf2)

            sigmaM = np.sqrt(sigma2_M)
            sigmaM_f = np.sqrt(sigma2_M_f)
            sigmaM_f2 = np.sqrt(sigma2_M_f2)

            results['D_z'][valid_mask] = Dz
            results['sigma_M'][valid_mask] = sigmaM
            results['sigma_M_f'][valid_mask] = sigmaM_f
            results['sigma_M_f2'][valid_mask] = sigmaM_f2

        return results if return_all else results['sigma_M']
