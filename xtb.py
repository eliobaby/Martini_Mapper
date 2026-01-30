from rdkit import Chem
from rdkit.Chem import AllChem
from xtb_md import run_xtb_md
import mdtraj as md
import numpy as np

def smiles_to_xyz(smiles, xyz_file="mol.xyz"):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()

    with open(xyz_file, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n\n")
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            pos = conf.GetAtomPosition(i)
            f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

def smiles_to_pdb(smiles, pdb_file="mol.pdb"):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    Chem.MolToPDBFile(mol, pdb_file)
    return pdb_file

def smiles_to_ref(name, smiles, *, drop_first_frac: float = 0.2, xtb_md_kwargs: dict = None):
    """
    Build AA reference trajectory for a molecule:
      - RDKit embed + UFF optimize -> name.xyz and name.pdb
      - xTB MD -> name_md/xtb.trj
      - Convert to ref_{name}.gro / ref_{name}.xtc for downstream CG fitting

    Parameters
    ----------
    drop_first_frac : float
        Fraction (0..1) of *valid* frames to drop from the beginning of the trajectory
        (simple equilibration discard). Default 0.2 (drop first 20%).
    xtb_md_kwargs : dict
        Optional kwargs passed to run_xtb_md() (e.g., time_ps, temp_k, step_fs, preopt).
    """
    if xtb_md_kwargs is None:
        xtb_md_kwargs = {}

    smiles_to_xyz(smiles, f"{name}.xyz")
    out = run_xtb_md(
        xyz_path=f"{name}.xyz",
        workdir=f"{name}_md",
        **xtb_md_kwargs
    )
    smiles_to_pdb(smiles, f"{name}.pdb")

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
    if drop_first_frac is None:
        drop_first_frac = 0.0
    drop_first_frac = float(drop_first_frac)

    drop_n = int(traj.n_frames * drop_first_frac)
    if drop_n >= traj.n_frames:
        drop_n = max(0, traj.n_frames - 1)

    if drop_n > 0:
        print(f"[smiles_to_ref] Dropping first {drop_n}/{traj.n_frames} frames ({drop_first_frac*100:.1f}%).")
        traj = traj[drop_n:]

    if traj.n_frames == 0:
        raise RuntimeError(
            "No frames left after dropping equilibration. "
            "Reduce drop_first_frac or increase MD length."
        )

    traj.save_gro(f"ref_{name}.gro")
    traj.save_xtc(f"ref_{name}.xtc")

