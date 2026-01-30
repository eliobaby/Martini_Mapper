from rdkit import Chem
from rdkit.Chem import AllChem
import math

def get_heavy_atom_coords_and_mol(smiles: str, randomSeed: int = 42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = randomSeed
    ok = AllChem.EmbedMolecule(mol, params)
    if ok != 0:
        raise ValueError("RDKit failed to embed a 3D conformer for this SMILES.")

    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    atom_data = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        atom_data.append({
            "index": idx,
            "symbol": atom.GetSymbol(),
            "x": pos.x,
            "y": pos.y,
            "z": pos.z,
        })

    return mol, atom_data


def gro_text_from_atom_data(atom_data, compound_name: str, box=None) -> str:
    if box is None:
        box = (10.0, 10.0, 10.0)  # nm

    res_name = "res"
    n_atoms = len(atom_data)

    lines = []
    lines.append(f"1{compound_name}")
    lines.append(f"{n_atoms}")

    for count, atom in enumerate(atom_data, start=1):
        atom_nr = count
        elem = atom["symbol"]
        atom_name = f"{elem}{count}"

        # Å -> nm
        x = atom["x"] / 10.0
        y = atom["y"] / 10.0
        z = atom["z"] / 10.0

        line = (
            f"{1:5d}{res_name:<5s}{atom_name:>5s}{atom_nr:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
        lines.append(line)

    bx, by, bz = box
    lines.append(f"{bx:10.5f}{by:10.5f}{bz:10.5f}")
    return "\n".join(lines) + "\n"


def itp_text_from_smiles_mol(mol, atom_data, compound_name: str) -> str:
    heavy_idx_to_itp = {a["index"]: i + 1 for i, a in enumerate(atom_data)}

    coord_nm = {}
    for a in atom_data:
        coord_nm[a["index"]] = (a["x"] / 10.0, a["y"] / 10.0, a["z"] / 10.0)

    lines = []
    lines.append("[ moleculetype ]")
    lines.append("; name        nrexcl")
    lines.append("res           1")
    lines.append("")

    lines.append("[ atoms ]")
    header_fmt = "; {:>5s} {:>5s} {:>5s} {:>5s} {:>5s} {:>5s} {:>10s} {:>10s}"
    lines.append(header_fmt.format("nr", "type", "resnr", "resid", "atom", "cgnr", "charge", "mass"))

    atom_fmt = "{:7d} {:5s} {:5d} {:5s} {:5s} {:5d} {:10.3f} {:10.3f}"

    element_mass = {
        "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974,
        "F": 18.998, "Cl": 35.45, "Br": 79.904, "I": 126.90,
        "B": 10.81, "Si": 28.085, "H": 1.008
    }

    for i, a in enumerate(atom_data, start=1):
        nr = i
        atype = a["symbol"][:5]
        resnr = 1
        resid = "res"
        elem = a["symbol"]
        atom_name = f"{elem}{i}"
        cgnr = nr
        charge = 0.0
        mass = float(element_mass.get(a["symbol"], 0.0))
        lines.append(atom_fmt.format(nr, atype, resnr, resid, atom_name, cgnr, charge, mass))

    lines.append("")
    lines.append("[ bonds ]")
    lines.append(";  ai    aj   funct      b0(nm)        k")
    bond_fmt = "{:5d} {:5d} {:5d} {:12.5f} {:8d}"

    k = 20000
    funct = 1
    emitted = set()

    for b in mol.GetBonds():
        a1 = b.GetBeginAtomIdx()
        a2 = b.GetEndAtomIdx()

        if a1 not in heavy_idx_to_itp or a2 not in heavy_idx_to_itp:
            continue

        ai = heavy_idx_to_itp[a1]
        aj = heavy_idx_to_itp[a2]
        key = (min(ai, aj), max(ai, aj))
        if key in emitted:
            continue
        emitted.add(key)

        x1, y1, z1 = coord_nm[a1]
        x2, y2, z2 = coord_nm[a2]
        dx, dy, dz = x1 - x2, y1 - y2, z1 - z2
        b0 = math.sqrt(dx * dx + dy * dy + dz * dz)

        lines.append(bond_fmt.format(ai, aj, funct, b0, k))

    lines.append("")
    return "\n".join(lines)