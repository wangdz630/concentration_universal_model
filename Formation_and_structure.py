"""
We calculate the halo concentration for various cosmological simulations using our universal concentration 
model, given by c_vir = A * (D_zf * D_z * sigma_M)**alpha + B. We first present the relations for D_zf derived 
from the simulations. Subsequently, we derive the concentration of a halo of mass Mat redshift z. In cases where 
D_zf is not directly required, the halo formation epoch can be obtained from a semi-analytical model or another 
fitting formula for the halo mass accretion history (MAH).
"""
import numpy as np
import sys
import os
# Define parameters for z_f calculation models
D_zf_MODELS = {
    "LCDM": {'alpha': -1.59,  'beta': 0.10, 'gamma': 2.29},
    "SCDM": {'alpha': -1.62,  'beta': 0.13, 'gamma': 2.33},
    "OCDM": {'alpha': -1.60,  'beta': 0.12, 'gamma': 2.34},
    "SF_n2": {'alpha': -1.60,  'beta': 0.14, 'gamma': 2.36},
    "LWDM": {'alpha': -1.60,  'beta': 0.13, 'gamma': 2.34},
    "wCDM": {'alpha': -1.60,  'beta': 0.13, 'gamma': 2.34},
}

def universal_model(D_z, D_zf, sigma_M, A=6.42, alpha=1.18, B=2.60):
    """Calculate c_vir"""
    return A * (D_zf * D_z * sigma_M)**alpha + B

def calculate_D_zf(D_z, sigma_M, D_zf_model='LCDM'):
    """Calculate formation redshift z_f"""
    params = D_zf_MODELS[D_zf_model]
    return D_z*(params['alpha'] * (sigma_M*D_z)**params['beta'] + params['gamma'])

def process_halo_data(input_file, cosmology_name, z=0.0, D_zf_model='LCDM'):
    """Process halo data and return results"""
    
    # Read data
    data = np.genfromtxt(input_file, skip_header=1)
    Mvir = data[:, 0]    # First column is Mvir
    D_z = data[:, 1]
    sigma_M = data[:, 2]  # Third column is sigma_M*D(z)
    
    # Initialize results array
    results = np.zeros((len(Mvir), 5))
    
    # Calculate for each halo
    for i in range(len(Mvir)):
        D_zf = calculate_D_zf(D_z[i], sigma_M[i], D_zf_model=D_zf_model)
        c_vir = universal_model(D_z[i], D_zf, sigma_M[i])
        results[i] = [Mvir[i], sigma_M[i], D_z[i], D_zf, c_vir]
    
    return Mvir, D_z, results

def main():
    # Parse command line arguments
    input_file = sys.argv[1]
    cosmology = sys.argv[2]
    redshift = float(sys.argv[3])
    D_zf_model = sys.argv[4]
    
    # Set output filename
    output_file = f"halo_properties2_{redshift}_{cosmology}.dat"
    
    # Process data
    Mvir, D_z, results = process_halo_data(input_file, cosmology, redshift, D_zf_model)
    
    # Print results in specified order
    print("=" * 60)
    print(f"Input file: {os.path.basename(input_file)}")
    print(f"Total halos processed: {len(Mvir)}")
    print(f"Mvir range: {Mvir.min():.2e} - {Mvir.max():.2e} M_sun/h")
    print(f"Cosmology: {cosmology}")
    print(f"Redshift: z = {redshift}")
    print(f"Linear growth factor D(z): {D_z[0]:.5f}")
    print(f"D(z_f,peak) model: {D_zf_model}")
    print(f"D(z_f,peak) range: {results[:,3].min():.3f} - {results[:,3].max():.3f}")
    print(f"c_vir range: {results[:,4].min():.2f} - {results[:,4].max():.2f}")
    print(f"Output file: {output_file}")
    print("=" * 60)
    
    # Save results
    np.savetxt(
        output_file, results,
        fmt='%.6e %.6f %.6f %.6f %.6f',
        header='Mvir[M_sun/h] sigma(M,0) D(z) D(z_f,peak) c_vir',
        comments='')

if __name__ == '__main__':
    main()
