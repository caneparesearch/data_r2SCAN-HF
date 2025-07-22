#!/usr/bin/env python3
"""
Generate VASP inputs (POSCAR, POTCAR, KPOINTS, INCAR) for a static r2SCAN calculation
using outputs from a prior 'r2SCAN_relax' run.

Usage:
    python prepare_r2scan_static.py

Requirements:
    - pymatgen
    - Pesudopotentials setup so that pymatgen can access them.
    - A completed relaxation in ./r2SCAN_relax/

Note:
    - Try to check the symmetry of the relaxed structure. If close to a higher 
    symmmety they try to symmetrize the structure for more effecient further hybrid 
    calculations.
"""

import os
import sys
import glob

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints, Incar

# --------------------------
# User-configurable settings
# --------------------------
RELAX_DIR = "r2SCAN_relax"     # Relaxation run directory (must exist)
STATIC_DIR = "r2SCAN"          # Output directory for static run
# K-point grid options:
KAPPA = 700                    # Choose sufficiently large grid, the same grid would be used for global hybrid calculations
FORCE_GAMMA = True             # Force Gamma-centered grid
# Parallelization
KPAR = 4
#NCORES_PER_NODE = 64
#NUMBER_OF_NODES = 2
NCORE = 4 # NCORE should be a factor of NCORES_PER_NODE
#total_mpi_ranks = NCORE_PER_NODE*NUMBER_OF_NODES
#npar = total_mpi_ranks//(NCORE*KPAR)

# --------------------------
# Sanity checks & setup
# --------------------------
if not os.path.isdir(RELAX_DIR):
    print(f"Error: '{RELAX_DIR}' not found. First run the r2SCAN_relax calculation.")
    sys.exit(1)

if os.path.isdir(STATIC_DIR):
    print(f"Error: '{STATIC_DIR}' already exists. Remove it before proceeding.")
    sys.exit(1)

os.makedirs(STATIC_DIR, exist_ok=False)

# Copy POTCAR from relaxation
os.system(f"cp {RELAX_DIR}/POTCAR {STATIC_DIR}/")

# --------------------------
# Locate last restart folder
# --------------------------
restart_dirs = sorted(glob.glob(os.path.join(RELAX_DIR, "restart*/")))
if not restart_dirs:
    print(f"Error: No restart*/ folders found under '{RELAX_DIR}'.")
    sys.exit(1)

last_restart = restart_dirs[-1]
print(f"Using structure from: {last_restart}")

# Load CONTCAR
contcar_path = os.path.join(last_restart, "CONTCAR")
structure = Structure.from_file(contcar_path)

# --------------------------
# Generate POSCAR
# --------------------------
#print("Space group:", structure.get_space_group_info())
print()
print("Checking symmetry. Symmetrize if VASP dosen't recognize higher symmetry")
for symprec in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1, 5e-1]:
    print("symprec:", symprec, structure.get_space_group_info(symprec=symprec, angle_tolerance=5))
print()
print("Formula:", structure.composition.reduced_formula, "\n")

# Copy CONTCAR to POSCAR for static run
poscar_out = os.path.join(STATIC_DIR, "POSCAR")
os.system(f"cp {contcar_path} {poscar_out}")
print(f"Written POSCAR -> {poscar_out}\n")

# --------------------------
# Generate KPOINTS
# --------------------------
kpoints = Kpoints.automatic_density(
    structure,
    KAPPA,
    force_gamma=FORCE_GAMMA
)

kpoints_out = os.path.join(STATIC_DIR, "KPOINTS")
kpoints.write_file(kpoints_out)
print(f"Written KPOINTS -> {kpoints_out}\n")

# --------------------------
# Generate INCAR
# --------------------------
# Start from relaxation INCAR, then tweak for static run
incar = Incar.from_file(os.path.join(RELAX_DIR, "INCAR"))

incar["ALGO"] = "Normal"    # switch to normal electronic steps
incar["IBRION"] = -1        # no ionic moves
incar["ISIF"] = 0           # no stress/relaxation
incar["LWAVE"] = True       # write WAVECAR for restart
incar["KPAR"] = KPAR
incar["NCORE"] = NCORE
incar["NSW"] = 0            # zero ionic steps

# Optional: remove symmetry settings if needed
incar.pop("ISYM", None)
incar.pop("SYMPREC", None)

incar_out = os.path.join(STATIC_DIR, "INCAR")
incar.write_file(incar_out)
print(f"Written INCAR -> {incar_out}\n")

# --------------------------
# Display snippets of inputs
# --------------------------
print()
os.system(f"cat {incar_out}")
print()
os.system(f"cat {kpoints_out}")
print()
os.system(f"head -n 10 {poscar_out}")
print()
os.system(f"grep TITEL {STATIC_DIR}/POTCAR")
