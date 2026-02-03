import os
import shutil

base_dir = "."

for fname in os.listdir(base_dir):
    if not fname.endswith(".smiles"):
        continue

    mol_name = fname.replace(".smiles", "")
    mol_dir = os.path.join(base_dir, mol_name)

    # create folder per molecule
    os.makedirs(mol_dir, exist_ok=True)

    src = os.path.join(base_dir, fname)
    dst = os.path.join(mol_dir, fname)

    # move file
    shutil.move(src, dst)
