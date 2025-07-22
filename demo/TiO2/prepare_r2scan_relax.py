#!/usr/bin/env python3
"""
Generate VASP inputs (POSCAR, POTCAR, KPOINTS, INCAR) for an r2SCAN relaxation.

Usage:
    python prepare_vasp_inputs.py

Requirements:
    - pymatgen
    - A CIF file in the current directory
"""

import os

from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Potcar
from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# --------------------------
# User‑configurable settings
# --------------------------
MAGMOM_STR = "6*0.6"                   # Magnetic moments string for INCAR, set it here or set it in the generated incar
OUT_DIR = "r2SCAN_relax"               # Output directory for VASP inputs
CIF_PATH = "./TiO2.cif"                # Input structure
FUNC = "PBE_54"                        # POTCAR functional
ISPIN = 1                              # Set to 1 is spin unpolarized, 2 if spin polarized

# Other settings to considered in the appropirate sections
# Set appropiate supercell matrix need to describe the magnetic ordering.
# Set the type of PAW potentials to be used sv or pv


# --------------------------
# Prepare output directory
# --------------------------
os.makedirs(OUT_DIR, exist_ok=False)

# --------------------------
# Load and preprocess structure
# --------------------------
# Read CIF, print space group and formula
structure = Structure.from_file(CIF_PATH)
sg=SpacegroupAnalyzer(structure)
print("Space group:", sg.get_space_group_symbol(), "(", sg.get_space_group_number(), ")")
print("Reduced Formula:", structure.composition.reduced_formula, "\n")

# Convert to primitive cell for more efficient relaxations
primitive = sg.find_primitive()
print("Primitive cell:")
print(primitive, "\n")

# (Optional) build a supercell from the primitive cell
# Here we use the identity matrix (no change), but you can adjust as needed
supercell_matrix = [[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]]
primitive.make_supercell(supercell_matrix)

# Write POSCAR
poscar_path = os.path.join(OUT_DIR, "POSCAR")
primitive.to(filename=poscar_path)
print(f"Written POSCAR -> {poscar_path}\n")

# --------------------------
# Generate POTCAR
# --------------------------
# Build POTCAR sequence in site order, appending "_sv" for TMs
potcar_seq = []
prev_elem = None
for site in primitive:
    elem = str(site.specie.element)
    if elem != prev_elem:
        tag = f"{elem}_sv" if Element(elem).is_transition_metal else elem
        potcar_seq.append(tag)
        prev_elem = elem

potcar = Potcar(potcar_seq, functional=FUNC)
potcar_path = os.path.join(OUT_DIR, "POTCAR")
potcar.write_file(potcar_path)
print(f"Written POTCAR ({FUNC}) -> {potcar_path}\n")

# --------------------------
# Generate KPOINTS
# --------------------------

KPOINT_LENGTH_DENSITY = 48             # k‑point length density for r2SCAN
FORCE_GAMMA = True                     # Force Gamma‑centered grid

print("Generating Gamma‑centered k‑point grid")
kpoints = Kpoints.automatic_density_by_lengths(
    primitive,
    [KPOINT_LENGTH_DENSITY] * 3,
    force_gamma=FORCE_GAMMA
)
kpoints_path = os.path.join(OUT_DIR, "KPOINTS")
kpoints.write_file(kpoints_path)
print(f"Written KPOINTS -> {kpoints_path}\n")

# --------------------------
# Generate INCAR
# --------------------------

KPAR = 16                               # k‑point parallelization
NCORE = 16                             # Core count per k‑point


incar_settings = {
    # Electronic & relaxation
    "ALGO": "All",
    "ENCUT": 700,
    "EDIFF": 1e-06,
    "SIGMA": 0.05,
    "PREC": "Accurate",
    "NELM": 200,
    "ISMEAR": 0,

    "EDIFFG": -0.01,
    "NSW": 99,
    "IBRION": 2,
    "ISIF": 3,

    # Spin & symmetry
    "ISPIN": ISPIN,
    "MAGMOM": MAGMOM_STR,
    "ISYM": 0,
    "SYMPREC": 1e-10,

    # Meta‑GGA & potentials
    "METAGGA": "R2scan",
    "LASPH": True,
    "LMAXMIX": 4,
    "LMIXTAU": True,

    # Parallelization
    "KPAR": KPAR,
    "NCORE": NCORE,

    # Outputs
    "LWAVE": False,
    "LCHARG": False,
    "LELF": False,
    "LVTOT": False,
    "LAECHG": False,

    # Misc
}

incar = Incar.from_dict(incar_settings)
incar_path = os.path.join(OUT_DIR, "INCAR")
incar.write_file(incar_path)
print(f"Written INCAR -> {incar_path}\n")

# --------------------------
# Display input snippets
# --------------------------
print()
os.system(f"cat {incar_path}")
print()
os.system(f"cat {kpoints_path}")
print()
os.system(f"head -n 10 {poscar_path}")
print()
os.system(f"grep TITEL {potcar_path}")
