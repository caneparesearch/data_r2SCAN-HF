#!/usr/bin/env python3
"""
Generate VASP inputs (POSCAR, POTCAR, KPOINTS, INCAR) for a series of hybrid-density
r2SCANX calculations, using WAVECAR (and other static inputs) from a prior r2SCAN run.

Usage:
    python prepare_r2scanx.py

Requirements:
    - pymatgen
    - Completed static r2SCAN run in ./r2SCAN/ (with WAVECAR, POTCAR, POSCAR, KPOINTS, INCAR)
    - ../submit_scripts/vasp_run_r2scan submission script
"""

import os
import sys

from pymatgen.io.vasp.inputs import Incar

# --------------------------
# User-configurable settings
# --------------------------
STATIC_DIR = "r2SCAN"                  # Directory of the completed r2SCAN static run
X_VALUES = [0.10, 0.25, 0.50, 0.75, 1.00, "HF"]  # Fraction of HF exchange (AEXX) for r2SCANX
KPAR = 4                               # k-point parallelization
NCORE = 4                             # cores per k-point

# --------------------------
# Sanity check
# --------------------------
if not os.path.isdir(STATIC_DIR):
    print(f"Error: '{STATIC_DIR}' not found. Run the r2SCAN static calculation first.")
    sys.exit(1)

# --------------------------
# Loop over hybrid fractions
# --------------------------
for X in X_VALUES:
    
    if X=="HF":
        hybrid_dir = "HF"
    else:
        pct = int(round(X * 100))
        hybrid_dir = f"{STATIC_DIR}{pct}"

    # Skip if already exists
    if os.path.isdir(hybrid_dir):
        print(f"SKIP: '{hybrid_dir}' already exists.")
        continue

    # Create directory for this hybrid-density run
    os.makedirs(hybrid_dir)
    print(f"\n=== Preparing r2SCANX (AEXX={X}) in '{hybrid_dir}' ===")

    # Copy static inputs and WAVECAR
    for fname in ("POTCAR", "POSCAR", "KPOINTS", "WAVECAR"):
        src = os.path.join(STATIC_DIR, fname)
        dst = os.path.join(hybrid_dir, fname)
        if os.path.exists(src):
            os.system(f"cp {src} {dst}")
        else:
            print(f"ERROR: {src} dosen't exist")
            sys.exit(1)
    print(f"Copied POTCAR, POSCAR, KPOINTS, WAVECAR -> {hybrid_dir}")

    # --------------------------
    # Generate INCAR
    # --------------------------
    incar = Incar.from_file(os.path.join(STATIC_DIR, "INCAR"))

    # Hybrid-density settings
    incar["LHFCALC"] = True       # turn on Hartree–Fock
    incar["AEXX"] = 1.0 if X=="HF" else X             # fraction of HF exchange
    if X == "HF":
        incar["AMGGAC"] = 0.0         # fraction of r2SCAN correlation
    else:
        incar["AMGGAC"] = 1.0         # fraction of r2SCAN correlation
    incar["PRECFOCK"] = "Fast"    # accelerate Fock build
    incar["LFOCKACE"] = True      # ACE acceleration for Fock
    incar["ALGO"] = "Normal"      # Davidson algorithm (recommended for hybrids)

    # Parallelization tweaks
    incar["KPAR"] = KPAR
    incar["NCORE"] = NCORE

    # Write INCAR
    incar_out = os.path.join(hybrid_dir, "INCAR")
    incar.write_file(incar_out)
    print(f"Written INCAR -> {incar_out}")

    # --------------------------
    # Quick sanity-check output
    # --------------------------
    os.system(f"grep -E 'LHFCALC|AEXX|AMGGAC|PRECFOCK|LFOCKACE' {incar_out}")
    os.system(f"head -n 5 {os.path.join(hybrid_dir, 'POSCAR')}")

    os.system(f"cp {STATIC_DIR}/vasp_run_r2x {hybrid_dir}/vasp_run")
    os.system(f"cd {hybrid_dir}/; sbatch vasp_run; cd ../")
print()
