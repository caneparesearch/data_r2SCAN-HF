# Demo: r<sup>2</sup>SCAN Relaxation, Static, Hybrid, and Non‑SCF r<sup>2</sup>SCANY\@r<sup>2</sup>SCANX Calculations for TiO<sub>2</sub>

This demonstration outlines a four-stage VASP workflow for TiO<sub>2</sub>, comprising:

1. **r<sup>2</sup>SCAN structural relaxation**
2. **r<sup>2</sup>SCAN static orbital generation**
3. **r<sup>2</sup>SCANX hybrid orbital calculations**
4. **Non-self-consistent r<sup>2</sup>SCANY\@r<sup>2</sup>SCANX energy evaluations**

Each stage is automated via a dedicated Python script that generates consistent VASP inputs, promoting reproducibility and enabling systematic parameter studies.

---

## Prerequisites

- **Software**
  - Python 3.8+ with:
    - `pymatgen` (≥ 2024.1)
    - `numpy`
  - VASP 6.4.2 or later, compiled with **r²SCAN** and hybrid-functional support
- **Data**
  - CIF file of a TiO<sub>2</sub> (e.g. `TiO2.cif`)
  - Magnetic ordering details (if applicable) for `MAGMOM` in `prepare_r2scan_relax.py`
  - PAW potentials available to `pymatgen` (e.g. `PBE_54`)
- **Scripts & Job Submission**
  - Python utilities:
    - `prepare_r2scan_relax.py`
    - `prepare_r2scan_static.py`
    - `prepare_r2scanx_static.py`
    - `prepare_r2scany@r2scanx.py`

---

## Workflow Breakdown

### 1. r<sup>2</sup>SCAN Structural Relaxation

**Script:** `prepare_r2scan_relax.py`

- Reads the CIF file; analyzes and symmetrizes the structure via `SpacegroupAnalyzer`
- Constructs primitive (and optional supercell) cells for magnetic ordering
- Generated inputs in `r2SCAN_relax/`:
  - `POSCAR` (from primitive/supercell)
  - `POTCAR` (TM potentials with `_sv` suffix)
  - Gamma-centered `KPOINTS` (length density = 48 $\AA^{-1}$, high density)
  - `INCAR` with:
    - `METAGGA = R2SCAN`
    - `ENCUT = 700 eV`, `EDIFF = 1e-6`, `EDIFFG = -0.01` (tight convergence)
    - `ISIF = 3` (full-relaxation)
    - `ISPIN`, `MAGMOM` (from experimental ordering)

```bash
python prepare_r2scan_relax.py
```

### 2. r<sup>2</sup>SCAN Static Orbital Generation

**Script:** `prepare_r2scan_static.py`

- Copies `POTCAR` from `r2SCAN_relax/`
- Retrieves final `CONTCAR` from `r2SCAN_relax/restart*/`
- Generated input in `r2SCAN/`:
  - `POSCAR` (from `CONTCAR`)
  - Gamma-centered `KPOINTS` (uniform density of 700 $\AA^{-3}$, low density)
  - `INCAR` with:
    - `ALGO = Normal`, `IBRION = -1`, `ISIF = 0`, `NSW = 0`
    - `LWAVE = .TRUE.` to generate `WAVECAR`
    - Meta-GGA settings carried over from relaxation

```bash
python prepare_r2scan_static.py
```

### 3. r<sup>2</sup>SCANX Hybrid-Density Calculations

**Script:** `prepare_r2scanx_static.py`

- Loops over X ∈ {0.10, 0.25, 0.50, 0.75, 1.00, HF}
- For each X:
  - Creates `r2SCAN{X*100}/` (or `HF/`)
  - Copies `POTCAR`, `POSCAR`, `KPOINTS`, `WAVECAR` from `r2SCAN/`
  - Updates `INCAR`:
    - `LHFCALC = .TRUE.`, `AEXX = X`, `AMGGAC = 1.0 (or 0.0 if X==HF)`
    - `PRECFOCK = Fast`, `LFOCKACE = .TRUE.`, `ALGO = Normal`

```bash
python prepare_r2scanx_static.py
```

### 4. Non-SCF r<sup>2</sup>SCANY\@r<sup>2</sup>SCANX Energy Evaluations

**Script:** `prepare_r2scany@r2scanx.py`

- Defines density X ∈ {0.10, 0.25, 0.50, 0.75, 1.00, HF} and functional Y ∈ {0.00, 0.05, 0.10, …, 0.95, 1.00, HF}
- Creates subdirectories under `r2SCANY@r2SCANXs/` named `r2SCAN{Y*100}@r2SCAN{X*100}/`
- Copies static inputs and `WAVECAR`
- Updates `INCAR`:
  - `LHFCALC = .TRUE.`, `AEXX = Y`, `AMGGAC = 1.0 (or 0.0 if Y==HF)`
  - `ALGO = Eigenval`, `NELM = 1`
  - `PRECFOCK = Fast`, `LFOCKACE = .TRUE.`
  - `LWAVE = .FALSE.`

```bash
python prepare_r2scany@r2scanx.py
```

---

## Directory Structure

```
├── prepare_r2scan_relax.py       # Relaxation input generator
├── prepare_r2scan_static.py      # Static input generator
├── prepare_r2scanx_static.py     # Hybrid‑density r²SCANX input generator
├── prepare_r2scany@r2scanx.py    # Non‑SCF r²SCANY@r²SCANX inputs generator
├── r2SCAN_relax/                 # r²SCAN relaxation directory
├── r2SCAN/                       # r²SCAN static run directory
├── r2SCAN{10,25,50,75,100}/      # Hybrid‑density directories
└── r2SCANY@r2SCANXs/             # Non‑SCF energy evaluation runs
    ├── r2SCANY5@r2SCAN25/          # Example: Y=0.05 @ X=0.25
    └── …
```

---