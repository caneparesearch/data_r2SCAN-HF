#!/usr/bin/env python3
"""
Generate VASP inputs (POSCAR, POTCAR, KPOINTS, WAVECAR, INCAR) for a series of non-self-consistent
r2SCANY@r2SCANX calculations at hybrid-density r2SCANX with r2SCANY, using WAVECAR from prior runs.

Usage:
    python prepare_r2scanx.py

Requirements:
    - pymatgen
    - A completed static r2SCAN calculation in ./r2SCAN/ (with WAVECAR, POTCAR, POSCAR, KPOINTS, INCAR)
    - A run script named "vasp_run" inside the STATIC_DIR
"""

import os
import sys
import shutil
from pymatgen.io.vasp.inputs import Incar

# --------------------------
# User-configurable settings
# --------------------------
STATIC_DIR = "r2SCAN"                     # Completed r2SCAN static-run dir
PARENT_DIR = "r2SCANY@r2SCANXs"            # Output parent directory
X_VALUES = [0.0, 0.10, 0.25, 0.50, 0.75, 1.00, "HF"]
Y_VALUES = [round(i/100, 2) for i in range(0, 101, 5)] + ["HF"]
KPAR = 4                                   # k-point parallelization
NCORE = 4                                  # cores per k-point

# --------------------------
# Sanity checks
# --------------------------
if not os.path.isdir(STATIC_DIR):
    print(f"Error: Static directory '{STATIC_DIR}' not found. Run the r2SCAN static calculation first.")
    sys.exit(1)

os.makedirs(PARENT_DIR, exist_ok=True)

# --------------------------
# Loop over hybrid-density and hybrid-functional fractions
# --------------------------
for X in X_VALUES:
    # Determine density directory
    if X == 0:
        density_dir = "r2SCAN"
    elif X == "HF":
        density_dir = "HF"
    else:
        pct_X = int(round(float(X) * 100,0))
        density_dir = f"r2SCAN{pct_X}"

    if not os.path.isdir(density_dir):
        print(f"SKIP: Density directory '{density_dir}' not found.")
        continue

    for Y in Y_VALUES:
        # Determine functional directory
        if Y == 0:
            func_dir = "r2SCAN"
        elif Y == "HF":
            func_dir = "HF"
        else:
            pct_Y = int(round(float(Y) * 100,0))
            func_dir = f"r2SCAN{pct_Y}"

        out_dir = os.path.join(PARENT_DIR, f"{func_dir}@{density_dir}")
        if os.path.isdir(out_dir):
            print(f"SKIP: '{out_dir}' already exists.")
            continue

        # Create output directory
        os.makedirs(out_dir)
        print(f"\n=== Preparing non-SCF run: {func_dir}@{density_dir} -> '{out_dir}' ===")

        # Copy static inputs
        for fname in ("POTCAR", "POSCAR", "KPOINTS", "WAVECAR"):
            src = os.path.join(density_dir, fname)
            dst = os.path.join(out_dir, fname)
            if os.path.exists(src):
                shutil.copy2(src, dst)
            else:
                print(f"ERROR: Source file '{src}' not found.")
                sys.exit(1)

        # Generate and write INCAR
        incar = Incar.from_file(os.path.join(density_dir, "INCAR"))
        incar["LHFCALC"] = True
        incar["AEXX"] = 1.00 if Y == "HF" else Y
        incar["AMGGAC"] = 0.0 if X == "HF" else 1.0
        incar["PRECFOCK"] = "Fast"
        incar["LFOCKACE"] = True
        incar["ALGO"] = "Eigenval"
        incar["NELM"] = 1
        incar["KPAR"] = KPAR
        incar["LWAVE"] = False
        incar["NCORE"] = NCORE
        incar.write_file(os.path.join(out_dir, "INCAR"))
        print(f"Written INCAR -> {out_dir}/INCAR")

        # Copy run script
        run_script="./vasp_run_r2yr2x"
        if os.path.exists(run_script):
            shutil.copy2(run_script, f"{out_dir}/vasp_run")
            os.system(f"cd {out_dir}; sbatch vasp_run; cd ../../")
        #else:
        #    print(f"WARNING: Run script '{run_script}' not found; please copy it manually.")

print("\nAll done.")
