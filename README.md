# Martini Mapper

A Python-based framework for automatically generating Martini 3 coarse-grained models directly from SMILES strings.  
This tool is designed to transform the complex, manual process of creating coarse-grained topologies into a fast, systematic, and reproducible workflow.

---

## Overview

The Martini 3 force field is a powerful tool for large-scale molecular simulations, but its use depends on a crucial *mapping* step where atoms are grouped into coarse-grained beads.  
This has traditionally been a slow, manual process requiring significant chemical expertise.

This framework automates that entire process. It uses a sophisticated, rule-based algorithm to intelligently partition a molecule and assign the correct bead types, producing simulation-ready files without any manual intervention.

---

## Key Features

- **Fully Automated**: Go from a SMILES string to simulation-ready GROMACS files in a single step.  
- **Rule-Based and Extensible**: Uses a detailed, internal dictionary of chemical fragments that can be expanded to map new structures.  
- **Hierarchical Mapping**: Employs a "most-constrained-first" strategy, mapping rigid ring systems before flexible chains to ensure a robust and chemically sound model.  
- **Standardized Outputs**: Generates GROMACS-compatible `.gro` (coordinate) and `.itp` (topology) files.  

---

## Installation

Because this project depends on **xtb (Python API)** and `pip` installs for `xtb` can fail on Windows (often requires building from source), the **recommended** setup is using **Anaconda/conda**.

### Option A (Recommended): Anaconda / conda via `environment.yaml`

1) Install Anaconda (if you don’t have it):
[Install Anaconda](https://docs.anaconda.com/free/anaconda/install/)

2) Create the environment and activate it:

```bash
# Clone the repository
git clone https://github.com/eliobaby/Martini_Mapper.git
cd Martini_Mapper

# Create the conda environment from the file
conda env create -f environment.yaml

# Activate
conda activate martini_mapper
```

3) Run:

```bash
cd martini_mapper
python main.py
```

> Notes:
> - `environment.yaml` uses **conda-forge** to install `xtb-python`, `rdkit`, and other compiled dependencies reliably.
> - If you update dependencies, update `environment.yaml` accordingly.

---

### Option B (Advanced): pip + venv (may fail for xtb on Windows)

> **Windows note:** If installation fails with Meson/MSVC errors like “Unknown compiler(s)” or “cl.exe not found”, `xtb` is being built from source and you need the MSVC C++ build tools.

```bash
# Clone the repository
git clone https://github.com/eliobaby/Martini_Mapper.git
cd Martini_Mapper

# Create a virtual environment
python -m venv .venv

# Activate it
# macOS/Linux:
source .venv/bin/activate
# Windows (PowerShell):
.venv\Scripts\Activate.ps1

# Upgrade pip
python -m pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

```

---

### Install the project (editable)

```bash
python -m pip install -e .

```

---

## Usage

### Run as a module (recommended)

Non-interactive (pass `NAME` and `SMILES`):

```bash
python -m martini_mapper Benzene "c1ccccc1"
```

Interactive (no positional args → prompts for Name/Smiles):

```bash
python -m martini_mapper
```

---

## CLI arguments

Positional:

- `name` (optional): compound name (e.g., `Benzene`). If omitted, you will be prompted.
- `smiles` (optional): SMILES string (e.g., `c1ccccc1`). If omitted, you will be prompted.

Options:

- `--no-xtb`  
  Skip xTB coordinate/CG generation.

- `--no-files`  
  Do not write output files (`.txt/.gro/.itp`). Useful for “library mode” and tests.

- `--out-dir PATH`  
  Directory to place final outputs.  
  Default: `./<name>/`

- `--keep-intermediates`  
  Do not delete intermediate xTB/MD files.

Examples:

```bash
# Run without xTB and write no files (fast / test-friendly)
python -m martini_mapper Benzene c1ccccc1 --no-xtb --no-files

# Put outputs in a specific folder
python -m martini_mapper Benzene c1ccccc1 --out-dir outputs/Benzene

# Keep intermediate xTB/MD files for debugging
python -m martini_mapper Benzene c1ccccc1 --keep-intermediates
```

---

## Testing

Install pytest and run the test suite from the repository root:

```bash
python -m pip install pytest
pytest -q
```

More verbose output:

```bash
pytest -ra -vv
```

---

## Project layout (important)

- `martini_mapper/` — importable package code
- `tests/` — pytest test suite
- Run via `python -m martini_mapper` (preferred)

---

## Troubleshooting

### `ModuleNotFoundError` for sibling modules
When running as a package, imports must be **relative** inside `martini_mapper/`, e.g.:

```python
from .setup_mapping import ...
```

not:

```python
from setup_mapping import ...
```

### xTB issues
If xTB isn’t available on your machine or you want a quick run, use:

```bash
python -m martini_mapper NAME SMILES --no-xtb
```

---

## Output Files

After running, by default, the script will generate three files in the same directory:

- `<MOLECULE_NAME>.gro`: A GROMACS coordinate file containing the 3D positions of each coarse-grained bead.  
- `<MOLECULE_NAME>.itp`: A GROMACS topology file that defines the bead types, bonds, and other parameters for the simulation.  
- `<MOLECULE_NAME>.txt`: A detailed summary file showing how each atom was mapped, which bead it belongs to, and how the beads are connected. This file is very useful for debugging and validating the mapping.

---

## Limitations

This framework is under active development. The current version has the following limitations:

- The internal mapping dictionary has been primarily developed using molecules containing **Carbon, Oxygen, and Nitrogen**.  
- Support for other elements like **Sulfur, Phosphorus, and halogens** is still experimental.  
- The algorithm may fail on molecules containing highly unusual or complex chemical fragments that are not yet present in the dictionary.

---

## Future Directions

- **Expand the Dictionary**: Systematically test the framework on more diverse chemical datasets to add new rules and increase the mapping success rate.  
- **Integrate Machine Learning**: Develop a machine learning model to assist with the tuning of bead parameters, which is currently a manual process.  

---

## Citation

If you use this work in your research, please cite our upcoming paper:  

Bigting, K. V.; Nag, S.; An, Y. *Martini Mapper: An Automated Fragment-Based Framework for Developing Coarse-Grained Models within the Martini 3 Framework.* 

---

## License
