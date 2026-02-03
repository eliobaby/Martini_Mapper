from rdkit import Chem
from rdkit.Chem import AllChem

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
smiles = input("SMILES: ")
coords2 = get_coordinates_from_smiles2(smiles)
coords = get_coordinates_from_smiles(smiles)
print(coords)

plot_molecule(coords)
plot_molecule(coords2)
'''