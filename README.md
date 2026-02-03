# Martini Mapper - An Automated Martini 3 Coarse-Graining Framework

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

## Installation (pip)

### 1) Clone the repository

```bash
git clone https://github.com/eliobaby/Martini_Mapper.git
cd Martini_Mapper

### 2) Create and activate a virtual environment

**macOS/Linux**
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip

**Windows (PowerShell)**
```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip

### 3) Install dependencies

```bash
pip install -r requirements.txt
```

### 4) Install the project (editable)

```bash
pip install -e .
```

### Dependencies

We pin dependency versions in `requirements.txt` to improve reproducibility for benchmarks and the paper.

- RDKit
- NumPy
- MDTraj
- xtb

**Note on xtb:** If your workflow uses the `xtb` executable (command-line program), you must also install the xtb binaries and ensure `xtb` is available on your PATH.

## Usage

The script is run interactively from the command line.

```bash
python main.py
```

You will be prompted for:  
- **Name**: The name of your molecule.  
- **SMILES**: The SMILES string representation of your molecule.  

### Example Session

```bash
$ python main.py
Name: Aspirin
SMILES: CC(=O)OC1=CC=CC=C1C(=O)O
```

---

## Output Files

After running, the script will generate three files in the same directory:

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



