# Demo: r²SCAN Relaxation, Static, Hybrid, and Non‑SCF r²SCANY\@r²SCANX Calculations for TiO₂

This document provides a comprehensive demonstration of the multi‑stage VASP workflow applied to TiO₂, encompassing:

1. **Structural relaxation** using the r²SCAN functional
2. **Static orbital generation** with r²SCAN
3. **Hybrid‑density r²SCANX** calculations (AEXX = 0.10, 0.25, 0.50, 0.75, 1.00 and HF)
4. **Non‑self‑consistent (non‑SCF) r²SCANY\@r²SCANX** energy evaluations

Each stage relies on a dedicated Python utility to generate standardized VASP inputs, facilitating reproducibility and parameter control.

---

## Prerequisites

Before proceeding, ensure the following software and data are available:

- **Python 3.8+** with the following packages installed:
  - `pymatgen` (≥ 2024.1)

- **VASP (≥6.4.2)** compiled with r2SCAN
- **CIF file** of the target TiO₂ polymorph in the working directory (e.g. `TiO2.rutile.cif`)
- **Magnetic ordering information** for TiO₂ (if applicable), to set appropriate `MAGMOM` strings
- **PAW potentials** accesible to `pymatgen` (e.g. PBE\_54 PAW potentials)

---

## Workflow and Scripts

### 1. Structural Relaxation (r²SCAN)

**Script:** `prepare_r2scan_relax.py`

- Reads the input CIF to generate a primitive cell and (optionally) supercell for magnetic ordering
- Generates in `r2SCAN_relax/`:
  - `POSCAR` from the (super)cell
  - `POTCAR` with PBE\_54 PAW potentials (for transition‑metal `_sv` potentials are used)
  - Gamma‑centered `KPOINTS` (length density = 48/\AA^{-1})
  - `INCAR` configured with:
    - `METAGGA = R2SCAN`, `LASPH = .TRUE.`
    - High cutoff (`ENCUT = 700 eV`), tight convergence (`EDIFF = 1e-6`)
    - Ionic relaxation parameters (`ISIF = 3`, `NSW = 99`, `IBRION = 2`)
    - Spin settings (`ISPIN`, `MAGMOM`)
    - Parallelization (`KPAR`, `NCORE`)

```bash
python prepare_r2scan_relax.py  # output: r2SCAN_relax/
```

### 2. Static r²SCAN (Orbital Generation)

**Script:** `prepare_r2scan_static.py`

- Copies `POTCAR` from `r2SCAN_relax/`
- Extracts the last `CONTCAR` from `r2SCAN_relax/restart*/`
- Writes in `r2SCAN/`:
  - `POSCAR` (from `CONTCAR`)
  - `KPOINTS` at density κ = 700/\AA^{-3}(Gamma‑centered)
  - Modified `INCAR`:
    - `ALGO = Normal`, `IBRION = -1`, `ISIF = 0`, `NSW = 0`
    - `LWAVE = .TRUE.` to generate `WAVECAR`
    - Retains high‑quality DFT settings from relaxation

```bash
python prepare_r2scan_static.py  # output: r2SCAN/
```

### 3. Hybrid‑Density r²SCANX Calculations

**Script:** `prepare_r2scanx_static.py`

- Iterates over HF‑exchange fractions: 0.10, 0.25, 0.50, 0.75, 1.00, and pure HF
- For each fraction X:
  - Creates `r2SCANX/` (e.g. `r2SCAN25/`)
  - Copies `POTCAR`, `POSCAR`, `KPOINTS`, and `WAVECAR` from `r2SCAN/`
  - Generates `INCAR` with hybrid settings:
    - `LHFCALC = .TRUE.`, `AEXX = X`, `AMGGAC = (X==HF ? 0.0 : 1.0)`
    - `PRECFOCK = Fast`, `LFOCKACE = .TRUE.`, `ALGO = Normal`
    - r2SCAN orbitals are used as the initial guess to speed up convergence.
    - Parallelization (`KPAR`, `NCORE`)

```bash
python prepare_r2scanx_static.py  # output: r2SCAN10/, r2SCAN25/, …, HF/
```

### 4. Non‑SCF r²SCANY\@r²SCANX Energy Evaluations

**Script:** `prepare_r2scany@r2scanx.py`

- Defines two parameter sets:
  - **Density** X ∈ {0.0, 0.10, …, 1.00, HF}
  - **Functional** Y ∈ {0.00, 0.05, …, 1.00, HF}
- Creates subfolders under `r2SCANY@r2SCANX/` named `{func_dir}@{density_dir}`
- Copies static inputs and `WAVECAR` from each density directory
- Writes `INCAR` with non‑SCF settings to evaluate the eigenvalues:
  - `LHFCALC = .TRUE.`, `AEXX = Y`, `AMGGAC = (X==HF ? 0.0 : 1.0)`
  - `ALGO = Eigenval`, `NELM = 1`, `PRECFOCK = Fast`, `LFOCKACE = .TRUE.`
  - Disables `LWAVE`, retains `KPAR`/`NCORE`

```bash
python "prepare_r2scany@r2scanx.py"  # output: r2SCANY@r2SCANX/
```

---

## Directory Structure

```
├── prepare_r2scan_relax.py       # Relaxation input generator
├── prepare_r2scan_static.py      # Static input generator
├── prepare_r2scanx_static.py     # Hybrid‑density r²SCANX inputs
├── prepare_r2scany@r2scanx.py    # Non‑SCF r²SCANY@r²SCANX inputs
├── r2SCAN_relax/                 # r²SCAN relaxation outputs
├── r2SCAN/                       # r²SCAN static run (WAVECAR, INCAR, etc.)
├── r2SCAN{10,25,50,75,100}/       # Hybrid‑density directories
└── r2SCANY@r2SCANXs/              # Non‑SCF energy evaluation runs
    ├── r2SCANY5@r2SCAN25/       # Example: Y=0.05 @ X=0.25
    └── …
```

---


