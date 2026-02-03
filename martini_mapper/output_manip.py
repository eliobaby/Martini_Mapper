from rdkit import Chem
from rdkit.Chem import AllChem
from coord_plotter import plot_molecule
def get_coordinates_from_smiles(smiles: str):
    """
    Given a SMILES string, generate a 3D conformer and return
    a list of atom coordinates in the order of the atom indices.

    Returns:
        List[Dict]: each dict contains {'index', 'symbol', 'x', 'y', 'z'}
    """
    # Parse SMILES and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol)

    # Embed 3D coordinates and optimize geometry
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    coords = []
    # Iterate atoms in RDKit atom order (matches SMILES order for heavy atoms)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        coords.append({
            'index': idx,
            'symbol': atom.GetSymbol(),
            'x': pos.x,
            'y': pos.y,
            'z': pos.z
        })
    return coords
'''
from openbabel import pybel

def get_coordinates_from_smiles2(smiles: str):
    """
    Given a SMILES string, generate a 3D conformer using Open Babel and return
    a list of heavy-atom coordinates in the order of the atom indices,
    ignoring all hydrogens.

    Returns:
        List[Dict]: each dict contains {'index', 'symbol', 'x', 'y', 'z'}
    """
    # Read SMILES and add hydrogens
    mol = pybel.readstring("smi", smiles)
    mol.addh()

    # Generate 3D coordinates and optimize geometry
    mol.make3D()
    ff = pybel._forcefields["uff"]  # UFF force field
    if not ff.Setup(mol.OBMol):
        raise RuntimeError("Failed to set up force field for optimization.")
    ff.ConjugateGradients(500, 1.0e-4)
    ff.GetCoordinates(mol.OBMol)

    coords = []
    # Iterate over atoms, skipping hydrogens
    for atom in mol.atoms:
        if atom.atomicnum == 1:
            continue
        # OBAtom idx is 1-based; subtract 1 for zero-based indexing
        idx = atom.idx - 1
        x, y, z = atom.coords
        # Extract element symbol from OBAtom type (letters before digits/dots)
        raw_type = atom.OBAtom.GetType()
        # take leading alphabetic characters as symbol
        symbol = ''
        for c in raw_type:
            if c.isalpha():
                symbol += c
            else:
                break
        # ensure proper casing (first uppercase, rest lowercase)
        symbol = symbol[0].upper() + symbol[1:].lower() if len(symbol) > 1 else symbol.upper()
        coords.append({
            'index': idx,
            'symbol': symbol,
            'x': x,
            'y': y,
            'z': z
        })
    return coords

smiles = input("SMILES: ")
coords2 = get_coordinates_from_smiles2(smiles)
coords = get_coordinates_from_smiles(smiles)
print(coords)

plot_molecule(coords)
plot_molecule(coords2)
'''