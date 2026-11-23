"""
Halo Analysis Pipeline
======================

Unified interface for the halo analysis pipeline combining:
  - Paper I   (Statistical halo concentration: peak, mean, median -> ./output_halo_properties/)
  - Paper II  (Individual halo c & concentration distributions -> ./output_individual_halo_properties/)
  - Paper III (Universal mass accretion histories, forward c(z) tracking,
               and Bayesian inverse reconstruction -> ./output_mah_properties/)

Usage:
  python halo_pipeline.py [cosmology] [redshift] [mode] [halo_def]

Arguments:
  cosmology : Cosmology model (planck18, WMAP7, etc.)
  redshift  : Redshift value (e.g., 0.0)
  mode      : Execution mode (optional, default: peak)
  halo_def  : Halo definition: vir, 200c, or both (optional, default: both)

Modes:
  - Paper I (Statistical concentrations):
      peak (default), mean, median, all
  - Paper II (Individual halo c & distributions):
      relaxed, all_halos, relaxed_einasto, individual
  - Paper III (Universal MAHs & structural evolution):
      mah          : Universal MAH trajectories M(z)/M0 for halos in input_masses.dat
      c_evolution  : Continuous forward concentration evolution c(z)
      reconstruct  : Bayesian MAP inversion from (M_obs, c_obs) in input_mass_concentration.dat
  - All Papers:
      all_papers   : Runs peak (Paper I), relaxed (Paper II), and mah (Paper III)

Examples:
  python halo_pipeline.py planck18 0.0                     # Paper I peak (both 200c and vir)
  python halo_pipeline.py planck18 0.0 peak vir            # Paper I peak for vir only
  python halo_pipeline.py planck18 0.0 relaxed 200c        # Paper II relaxed for 200c only
  python halo_pipeline.py planck18 0.0 all_halos vir       # Paper II all halos for vir only
  python halo_pipeline.py planck18 0.0 mah vir             # Paper III MAH for vir only
  python halo_pipeline.py planck18 0.0 c_evolution 200c    # Paper III c(z) for 200c only
  python halo_pipeline.py planck18 0.0 reconstruct vir     # Paper III MAP reconstruction for vir only
  python halo_pipeline.py planck18 0.0 mah                 # Paper III MAH (both 200c and vir)
  python halo_pipeline.py planck18 0.0 all_papers          # All papers (both 200c and vir)
"""
import sys
import subprocess
import os
import time
import warnings
from datetime import datetime
import numpy as np

def main():
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print(__doc__)
        sys.exit(1)
    
    # halo mass file
    mass_file = "input_masses.dat"
    reconstruct_file = "input_mass_concentration.dat"
    
    cosmology = sys.argv[1]
    redshift = sys.argv[2]
    
    mode_input = sys.argv[3].lower() if len(sys.argv) >= 4 else 'peak'
    
    paper1_modes = ['peak', 'mean', 'median']
    paper2_modes = ['relaxed', 'all_halos', 'relaxed_einasto']
    paper3_modes = ['mah', 'c_evolution', 'reconstruct']

    if mode_input == 'all':
        modes_to_run = paper1_modes
    elif mode_input == 'individual':
        modes_to_run = paper2_modes
    elif mode_input == 'all_papers':
        modes_to_run = ['peak', 'relaxed', 'mah', 'c_evolution', 'reconstruct']
    elif mode_input in paper1_modes or mode_input in paper2_modes or mode_input in paper3_modes:
        modes_to_run = [mode_input]
    else:
        print(
            f"ERROR: Invalid mode '{mode_input}'.\n"
            "Choose from:\n"
            "  Paper I  : peak, mean, median, all\n"
            "  Paper II : relaxed, all_halos, relaxed_einasto, individual\n"
            "  Paper III: mah, c_evolution, reconstruct\n"
            "  Combined : all_papers"
        )
        sys.exit(1)
        
    # Determine halo boundary definitions to evaluate
    if len(sys.argv) == 5:
        def_choice = sys.argv[4].lower()
        if def_choice in ['vir', '200c']:
            halo_defs = [def_choice]
        elif def_choice in ['both', 'all', 'all_defs']:
            halo_defs = ['200c', 'vir']
        else:
            print(f"ERROR: Invalid halo definition '{def_choice}'. Choose from: vir, 200c, or both.")
            sys.exit(1)
    else:
        # Default behavior: run both definitions
        halo_defs = ['200c', 'vir']
    
    # Start timing
    start_time = time.time()
    start_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    print("=" * 65)
    print("UNIFIED HALO ANALYSIS PIPELINE (PAPERS I, II, & III)")
    print("=" * 65)
    print(f"Start time: {start_datetime}")
    print(f"Input Data: {mass_file}")
    print(f"Cosmology:  {cosmology}")
    print(f"Redshift:   z = {redshift}")
    print(f"Modes:      {', '.join(modes_to_run)}")
    print("=" * 65)
    
    generated_files = []
    
    # Check if any Paper I or Paper II modes need to be run
    run_p1_p2 = any(m in paper1_modes or m in paper2_modes for m in modes_to_run)
    
    if run_p1_p2:
        # ---------------------------------------------------------
        # Step 1: Compute halo mass variance
        # ---------------------------------------------------------
        print("\n[STEP 1] Computing halo mass variance...")
        cmd1 = [sys.executable, "halo_variance.py", mass_file, cosmology, redshift]
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
        
        # ---------------------------------------------------------
        # Step 2: Compute halo concentrations (Paper I & Paper II)
        # ---------------------------------------------------------
        print("\n[STEP 2] Computing halo concentrations (Paper I & Paper II)...")
        for current_mode in modes_to_run:
            if current_mode in paper3_modes:
                continue
            if current_mode == 'relaxed_einasto':
                if '200c' not in halo_defs:
                    print(f"SKIPPING: relaxed_einasto is calibrated for 200c only (Paper II Eq. 8).")
                    continue
                defs_for_mode = ['200c']
            else:
                defs_for_mode = halo_defs
            for current_def in defs_for_mode:
                print(f"\n--- Processing: Mode = {current_mode.upper()}, Def = {current_def.upper()} ---")
                cmd2 = [sys.executable, "halo_model.py", "output_variance.dat", cosmology, redshift, current_def, current_mode]
                print(f"Running: {' '.join(cmd2)}")
                result2 = subprocess.run(cmd2)
                if result2.returncode != 0:
                    print(f"ERROR: Step 2 failed for {current_def} {current_mode}")
                    sys.exit(1)
                    
                if current_mode in paper2_modes:
                    expected_output = f"./output_individual_halo_properties/c_distribution_{float(redshift)}_{cosmology}_{current_def}_{current_mode}.dat"
                else:
                    expected_output = f"./output_halo_properties/halo_properties_{float(redshift)}_{cosmology}_{current_def}_{current_mode}.dat"

                if os.path.exists(expected_output):
                    generated_files.append(expected_output)
                else:
                    print(f"WARNING: Expected output {expected_output} not found.")

    # ---------------------------------------------------------
    # Step 3: Paper III Tasks (MAH, c_evolution, reconstruct)
    # ---------------------------------------------------------
    for current_mode in modes_to_run:
        if current_mode not in paper3_modes:
            continue
        print(f"\n[PAPER III] Executing {current_mode.upper()}...")
        for current_def in halo_defs:
            print(f"\n--- Processing Paper III: Task = {current_mode.upper()}, Def = {current_def.upper()} ---")
            if current_mode == 'mah':
                cmd3 = [sys.executable, "mah_model.py", "mah",
                        "--input", mass_file,
                        "--cosmology", cosmology,
                        "--z", str(redshift),
                        "--halo-def", current_def]
                print(f"Running: {' '.join(cmd3)}")
                result3 = subprocess.run(cmd3)
                if result3.returncode != 0:
                    print(f"ERROR: Paper III {current_mode} failed for {current_def}")
                    sys.exit(1)
            elif current_mode == 'c_evolution':
                for rel in ['relaxed', 'all_halos']:
                    cmd3 = [sys.executable, "mah_model.py", "c_evolution",
                            "--input", mass_file,
                            "--cosmology", cosmology,
                            "--z", str(redshift),
                            "--halo-def", current_def,
                            "--relation", rel]
                    print(f"Running: {' '.join(cmd3)}")
                    result3 = subprocess.run(cmd3)
                    if result3.returncode != 0:
                        print(f"ERROR: Paper III {current_mode} failed for {current_def} ({rel})")
                        sys.exit(1)
            elif current_mode == 'reconstruct':
                input_file = reconstruct_file if os.path.exists(reconstruct_file) else mass_file
                if not os.path.exists(input_file):
                    print(f"WARNING: Reconstruction input {input_file} not found, skipping.")
                    continue
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        chk = np.loadtxt(input_file, max_rows=2, comments='#', ndmin=2)
                    if chk.shape[1] < 2:
                        print(f"SKIPPING reconstruct: '{input_file}' needs at least 2 columns (M_obs c_obs).")
                        continue
                except Exception:
                    continue
                for rel in ['relaxed', 'all_halos']:
                    cmd3 = [sys.executable, "mah_model.py", "reconstruct",
                            "--input", input_file,
                            "--cosmology", cosmology,
                            "--z", str(redshift),
                            "--halo-def", current_def,
                            "--relation", rel]
                    print(f"Running: {' '.join(cmd3)}")
                    result3 = subprocess.run(cmd3)
                    if result3.returncode != 0:
                        print(f"ERROR: Paper III {current_mode} failed for {current_def} ({rel})")
                        sys.exit(1)

    # Collect newly generated output_mah_properties files
    if os.path.exists("./output_mah_properties"):
        for fname in sorted(os.listdir("./output_mah_properties")):
            if fname.endswith(".dat"):
                fpath = os.path.join("./output_mah_properties", fname)
                if os.path.getmtime(fpath) >= start_time - 1.0 and fpath not in generated_files:
                    generated_files.append(fpath)

    total_time = time.time() - start_time
    
    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------
    if generated_files:
        print("\n" + "=" * 65)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 65)
        print("Output files available:")
        for f in generated_files[:15]:
            if os.path.exists(f):
                file_size = os.path.getsize(f)
                print(f"  - {f} ({file_size:,} bytes)")
        if len(generated_files) > 15:
            print(f"  ... and {len(generated_files) - 15} more files.")
        print(f"\nTotal time elapsed: {total_time:.1f} seconds")
        print("=" * 65)
    else:
        print("Pipeline finished.")

if __name__ == "__main__":
    main()
