from rdkit import Chem
from rdkit.Chem import AllChem
from xtb_md_runner import run_xtb_md
import mdtraj as md
import numpy as np

def build_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    return mol

def write_xyz_from_mol(mol, xyz_file):
    conf = mol.GetConformer()
    with open(xyz_file, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n\n")
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            pos = conf.GetAtomPosition(i)
            f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

def write_pdb_from_mol(mol, pdb_file):
    Chem.MolToPDBFile(mol, pdb_file)

def smiles_to_ref(name, smiles, *, drop_first_frac: float = 0):
    """
    Build AA reference trajectory for a molecule:
      - RDKit embed + UFF optimize -> name.xyz and name.pdb
      - xTB MD -> name_md/xtb.trj
      - Convert to ref_{name}.gro / ref_{name}.xtc for downstream CG fitting

    Parameters
    ----------
    drop_first_frac : float
        Fraction (0..1) of *valid* frames to drop from the beginning of the trajectory
        (simple equilibration discard). Default 0.8 (drop first 80%).
    """
    mol = build_mol(smiles)

    xyz_path = f"{name}.xyz"
    pdb_path = f"{name}.pdb"

    write_xyz_from_mol(mol, xyz_path)
    write_pdb_from_mol(mol, pdb_path)

    out = run_xtb_md(
        xyz_path=xyz_path,
        workdir=f"{name}_md",
    )

    traj = md.load_xyz(f"{name}_md/xtb.trj", top=f"{name}.pdb")

    # ---- NaN/Inf guard + frame filter ----
    finite_frame = np.isfinite(traj.xyz).all(axis=(1,2))
    if not finite_frame.all():
        bad = int((~finite_frame).sum())
        raise ValueError(f"[smiles_to_ref] {bad} frames contain NaN/Inf. Aborting.")

    if traj.n_frames == 0:
        raise RuntimeError(
            "All frames were invalid (NaN/Inf). Check xtb logs:\n"
            f"  {out['stdout_log']}\n"
            f"  {out['stderr_log']}"
        )

    # ---- Drop equilibration frames (first X%) ----
    drop_n = int(traj.n_frames * float(drop_first_frac))
    drop_n = min(drop_n, traj.n_frames - 1)

    if drop_n > 0:
        traj = traj[drop_n:]

    traj.save_gro(f"ref_{name}.gro")
    traj.save_xtc(f"ref_{name}.xtc")
