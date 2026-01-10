"""
Halo Analysis Pipeline
======================

Unified interface for the two-step halo analysis pipeline.

Usage: python halo_pipeline.py [mass_file] [cosmology] [redshift]
Example: python halo_pipeline.py input_masses.dat LCDM 0.0
"""
import sys
import subprocess
import os
import time
from datetime import datetime

def main():
    if len(sys.argv) != 4:
        print(__doc__)
        print("\nArguments:")
        print("  mass_file  : Path to halo mass file")
        print("  cosmology  : Cosmology model (LCDM, etc.)")
        print("  redshift   : Redshift value")
        print("\nExample: python halo_pipeline.py input_masses.dat LCDM 0.0")
        sys.exit(1)
    
    mass_file, cosmology, redshift = sys.argv[1], sys.argv[2], sys.argv[3]
    
    # Start timing
    start_time = time.time()
    start_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    print("=" * 50)
    print("HALO ANALYSIS PIPELINE")
    print("=" * 50)
    print(f"Start time: {start_datetime}")
    print(f"Input:      {mass_file}")
    print(f"Cosmology:  {cosmology}")
    print(f"Redshift:   z = {redshift}")
    print("=" * 50)
    
    # Step 1: Compute halo mass variance
    print("\n[STEP 1/2] Computing halo mass variance...")
    cmd1 = ["python", "halo_variance.py", mass_file, cosmology, redshift]
    print(f"Running: {' '.join(cmd1)}")
    
    step1_start = time.time()
    result1 = subprocess.run(cmd1)
    step1_time = time.time() - step1_start
    
    if result1.returncode != 0:
        print("ERROR: Step 1 failed")
        sys.exit(1)
    
    if not os.path.exists("output_variance.dat"):
        print("ERROR: output_variance.dat not generated")
        sys.exit(1)
    
    print(f"Step 1 completed ({step1_time:.1f} seconds)")
    
    # Step 2: Compute halo concentrations
    print("\n[STEP 2/2] Computing halo concentrations...")
    cmd2 = ["python", "Formation_and_structure.py", "output_variance.dat", cosmology, redshift]
    print(f"Running: {' '.join(cmd2)}")
    
    step2_start = time.time()
    result2 = subprocess.run(cmd2)
    step2_time = time.time() - step2_start
    
    if result2.returncode != 0:
        print("ERROR: Step 2 failed")
        sys.exit(1)
    
    # Check output file
    output_files = [
        f"halo_properties_z{redshift}_{cosmology}.dat",
        f"halo_properties_{redshift}_{cosmology}.dat"
    ]
    
    output_file = None
    for f in output_files:
        if os.path.exists(f):
            output_file = f
            break
    
    total_time = time.time() - start_time
    
    if output_file:
        print("\n" + "=" * 50)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 50)
        print(f"Results saved to: {output_file}")
        file_size = os.path.getsize(output_file)
        print(f"File size:       {file_size:,} bytes")
        print("\nTIMING SUMMARY:")
        print(f"Step 1 (variance):    {step1_time:.1f} seconds")
        print(f"Step 2 (concentration): {step2_time:.1f} seconds")
        print(f"Total time:           {total_time:.1f} seconds")
        print("=" * 50)
    else:
        print("ERROR: Could not find output file")
        print(f"Total time elapsed: {total_time:.1f} seconds")
        sys.exit(1)

if __name__ == "__main__":
    main()