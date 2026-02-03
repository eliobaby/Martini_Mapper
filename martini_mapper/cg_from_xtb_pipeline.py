from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Set

import numpy as np
import mdtraj as md
from rdkit import Chem

import math

# -----------------------------
# Angle / dihedral helpers
# -----------------------------

kB = 0.008314462618  # kJ/mol/K
RAD2DEG = 180.0 / np.pi
DEG2RAD = np.pi / 180.0

def _wrap_to_pi(x: np.ndarray) -> np.ndarray:
    """Wrap radians to (-pi, pi]."""
    return (x + np.pi) % (2.0 * np.pi) - np.pi

def _circular_mean(phi: np.ndarray) -> float:
    """Circular mean for angles in radians."""
    return float(np.arctan2(np.sin(phi).mean(), np.cos(phi).mean()))

def _circular_std(phi: np.ndarray) -> float:
    """Circular standard deviation in radians: sqrt(-2 ln R)."""
    s = np.sin(phi).mean()
    c = np.cos(phi).mean()
    R = float(np.sqrt(s*s + c*c))
    R = max(R, 1e-300)  # numerical safety only
    return float(np.sqrt(max(0.0, -2.0 * np.log(R))))

# Reuse your existing functions from text_file.py
# (adjust import path/name as needed)
from text_file import (
    group_beads_by_type,
    compute_beads_connections,
    fix_beadtypes,
)

# -----------------------------
# 1) RDKit: heavy atom -> attached hydrogen indices
# -----------------------------

def build_heavy_to_hydrogens(smiles: str) -> Dict[int, List[int]]:
    """
    RDKit atom indexing for Chem.AddHs(mol):
      - Heavy atoms keep their original indices (0..nHeavy-1)
      - Explicit hydrogens are appended after heavy atoms

    Returns dict: heavy_index -> [H_index, H_index, ...]
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    heavy_to_h: Dict[int, List[int]] = {}

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        h_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]
        heavy_to_h[atom.GetIdx()] = h_neighbors

    return heavy_to_h


def expand_beads_with_hydrogens(
    bead_heavy_atoms: List[List[int]],
    heavy_to_h: Dict[int, List[int]],
) -> List[List[int]]:
    """
    For each bead group (heavy atoms only), include all attached H indices.
    """
    bead_atoms: List[List[int]] = []
    for heavy_list in bead_heavy_atoms:
        s: Set[int] = set(heavy_list)
        for hidx in heavy_list:
            for hid in heavy_to_h.get(hidx, []):
                s.add(hid)
        bead_atoms.append(sorted(s))
    return bead_atoms


# -----------------------------
# 2) Map AA trajectory -> bead trajectory (center of geometry)
# -----------------------------

def map_aa_to_cg_traj(
    ref_gro: str,
    ref_xtc: str,
    bead_atom_indices: List[List[int]],
    bead_names: List[str],
    out_cg_gro: str,
    out_cg_xtc: str,
) -> md.Trajectory:
    """
    Loads AA trajectory and writes CG trajectory where each bead position is
    mean(xyz) across its atom indices (center of geometry).
    """
    aa = md.load(ref_xtc, top=ref_gro)
    if not np.isfinite(aa.xyz).all():
        raise RuntimeError("AA trajectory contains NaN/Inf; ref.gro/ref.xtc are corrupted.")


    n_frames = aa.n_frames
    n_beads = len(bead_atom_indices)

    cg_xyz = np.zeros((n_frames, n_beads, 3), dtype=np.float32)
    for bi, idxs in enumerate(bead_atom_indices):
        idxs = np.asarray(idxs, dtype=int)
        cg_xyz[:, bi, :] = aa.xyz[:, idxs, :].mean(axis=1)

    # Minimal CG topology: 1 residue with n_beads atoms
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue("res", chain)
    for name in bead_names:
        top.add_atom(name[:4], element=md.element.carbon, residue=res)

    cg = md.Trajectory(cg_xyz, top)

    # carry unit cell if present
    if aa.unitcell_lengths is not None:
        cg.unitcell_lengths = aa.unitcell_lengths
        cg.unitcell_angles = aa.unitcell_angles

    cg.save_gro(out_cg_gro)
    cg.save_xtc(out_cg_xtc)
    return cg

def map_aa_to_cg_traj_from_loaded_aa(
    aa: md.Trajectory,
    bead_atom_indices: List[List[int]],
    bead_names: List[str],
    out_cg_gro: str,
    out_cg_xtc: str,
) -> md.Trajectory:
    n_frames = aa.n_frames
    n_beads = len(bead_atom_indices)

    cg_xyz = np.zeros((n_frames, n_beads, 3), dtype=np.float32)
    for bi, idxs in enumerate(bead_atom_indices):
        idxs = np.asarray(idxs, dtype=int)
        cg_xyz[:, bi, :] = aa.xyz[:, idxs, :].mean(axis=1)  # center of geometry

    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue("res", chain)
    for name in bead_names:
        top.add_atom(name[:4], element=md.element.carbon, residue=res)

    cg = md.Trajectory(cg_xyz, top)

    if aa.unitcell_lengths is not None:
        cg.unitcell_lengths = aa.unitcell_lengths
        cg.unitcell_angles = aa.unitcell_angles

    cg.save_gro(out_cg_gro)
    cg.save_xtc(out_cg_xtc)
    return cg

# -----------------------------
# 3) Fit bonds + angles from CG trajectory
# -----------------------------

def k_from_sigma(T: float, sigma: float) -> float:
    """Harmonic: k ≈ kBT / sigma^2. T in K. sigma in nm or radians."""
    kB = 0.008314462618  # kJ/mol/K
    if sigma <= 0:
        return 0.0
    return (kB * T) / (sigma * sigma)


def generate_angle_list_from_bonds(n_beads: int, bonds: List[Tuple[int, int]]) -> List[Tuple[int, int, int]]:
    """
    From bond graph, generate unique angles i-j-k (j is central).
    We'll emit each angle once using i<k to avoid duplicates.
    """
    neigh = {i: set() for i in range(n_beads)}
    for i, j in bonds:
        neigh[i].add(j)
        neigh[j].add(i)

    angles = set()
    for j in range(n_beads):
        nb = sorted(neigh[j])
        for a in range(len(nb)):
            for b in range(a + 1, len(nb)):
                i = nb[a]
                k = nb[b]
                # canonical: (min(i,k), j, max(i,k))
                angles.add((min(i, k), j, max(i, k)))

    return sorted(angles)


def fit_bonds_and_angles(
    cg: md.Trajectory,
    bonds: List[Tuple[int, int]],
    angles: List[Tuple[int, int, int]],
    T: float = 300.0,
):
    """
    Returns fitted parameters:
      bonds_fit: list of (i,j, b0_mean_nm, k_kJmol_nm2, std_nm)
      angles_fit: list of (i,j,k, th0_mean_rad, k_kJmol_rad2, std_rad)

    Notes
    -----
    - Bonds use k = kBT / sigma^2 (harmonic)
    - Angles use circular mean/std (more stable than linear stats)
    """
    bonds_fit = []
    for (i, j) in bonds:
        vals = md.compute_distances(cg, np.array([[i, j]], dtype=int))[:, 0]
        mu = float(vals.mean())
        #mu = float(vals[0]) #use this for first frame
        sigma = float(vals.std(ddof=1))
        k = float(k_from_sigma(T, sigma))
        bonds_fit.append((i, j, mu, k, sigma))

    angles_fit = []
    if len(angles) > 0:
        trip = np.array(angles, dtype=int)
        vals = md.compute_angles(cg, trip)  # radians
        for idx, (i, j, k_) in enumerate(angles):
            v = vals[:, idx].astype(float)
            th0 = _circular_mean(v)
            sigma = _circular_std(v)
            k = float(k_from_sigma(T, sigma))
            angles_fit.append((i, j, k_, th0, k, sigma))

    return bonds_fit, angles_fit

def generate_dihedral_list_from_bonds(n_beads: int, bonds: List[Tuple[int, int]]) -> List[Tuple[int,int,int,int]]:
    """
    Generate dihedrals i-j-k-l from bond graph:
      for each bond (j,k), choose neighbor i of j (i != k) and neighbor l of k (l != j).
    Deduplicate by canonical ordering of the tuple.
    """
    neigh = {i: set() for i in range(n_beads)}
    for a,b in bonds:
        neigh[a].add(b)
        neigh[b].add(a)

    dihs = set()
    for j,k in bonds:
        for i in neigh[j]:
            if i == k: 
                continue
            for l in neigh[k]:
                if l == j or l == i:
                    continue
                # canonicalize to avoid duplicates: pick lexicographically smallest between (i,j,k,l) and (l,k,j,i)
                t1 = (i,j,k,l)
                t2 = (l,k,j,i)
                dihs.add(min(t1,t2))
    return sorted(dihs)

def fit_dihedrals_simple(
    cg: md.Trajectory,
    dihedrals: List[Tuple[int,int,int,int]],
    T: float = 300.0,
):
    """
    Fit proper dihedrals from the mapped CG trajectory using Boltzmann inversion + Fourier fit.

    Produces parameters consistent with periodic dihedral form (GROMACS funct=1):
        U(φ) ≈ Σ_n K_n * (1 + cos(n*φ - φ0_n))

    Returns list of tuples: (i, j, k, l, phi0, K, n)
    Potentially multiple rows per dihedral (n = 1..N).
    """
    if not dihedrals:
        return []

    kBT = kB * T

    dih_arr = np.array(dihedrals, dtype=int)
    phis = md.compute_dihedrals(cg, dih_arr)  # (n_frames, n_dih), radians in [-pi, pi]

    nbins = 180
    n_terms = 3
    eps = 1e-12

    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Design matrix: [1, cos(φ), sin(φ), cos(2φ), sin(2φ), ...]
    X_cols = [np.ones_like(centers)]
    for n in range(1, n_terms + 1):
        X_cols.append(np.cos(n * centers))
        X_cols.append(np.sin(n * centers))
    X = np.vstack(X_cols).T  # (nbins, 1 + 2*n_terms)

    out: List[Tuple[int,int,int,int,float,float,int]] = []

    for col, (i, j, k, l) in enumerate(dihedrals):
        phi = _wrap_to_pi(phis[:, col].astype(float))

        hist, _ = np.histogram(phi, bins=edges, density=True)
        P = np.maximum(hist, eps)
        U = -kBT * np.log(P)
        U -= U.mean()  # remove constant offset

        coef, *_ = np.linalg.lstsq(X, U, rcond=None)

        for n in range(1, n_terms + 1):
            a_n = float(coef[1 + 2*(n-1)])
            b_n = float(coef[1 + 2*(n-1) + 1])

            K = math.sqrt(a_n*a_n + b_n*b_n)
            if K == 0.0:
                continue
            phi0 = math.atan2(b_n, a_n)
            out.append((i, j, k, l, float(phi0), float(K), int(n)))

    return out


# -----------------------------
# 4) Write .itp using your style (atoms + bonds + angles)
# -----------------------------

def bead_mass_from_type(bt: str) -> float:
    # same logic as export_bead_mapping(): T=36, S=54, else 72 :contentReference[oaicite:1]{index=1}
    if bt.startswith("T"):
        return 36.0
    if bt.startswith("S"):
        return 54.0
    return 72.0

def write_itp_with_bonds_angles_dihedrals(
    compound_name: str,
    bead_types_fixed: List[str],
    bonds_fit,
    angles_fit,
    dihedrals_fit,
    out_itp: str,
):
    """
    Writes a self-contained .itp for a single residue named 'res', similar to your itp_maker(),
    but uses fitted bond/angle values (instead of constant k=20000 and b0 from one frame).
    """
    lines = []
    lines.append("[ moleculetype ]")
    lines.append("; name        nrexcl")
    lines.append("res           1")
    lines.append("")

    # atoms
    lines.append("[ atoms ]")
    lines.append(";    nr  type  resnr resid  atom  cgnr     charge       mass")
    for idx, bt in enumerate(bead_types_fixed, start=1):
        mass = bead_mass_from_type(bt)
        lines.append(f"{idx:7d} {bt[:5]:5s} {1:5d} {'res':5s} {'C'+str(idx):5s} {idx:5d} {0.0:10.3f} {mass:10.3f}")
    lines.append("")

    # bonds (funct=1 harmonic)
    lines.append("[ bonds ]")
    lines.append(";  ai    aj   funct      b0(nm)          k(kJ/mol/nm^2)")
    for (i, j, b0, k, sigma) in bonds_fit:
        lines.append(f"{i+1:5d} {j+1:5d} {1:5d} {b0:12.5f} {k:16.1f}")
    lines.append("")

    # angles (funct=1 harmonic, in radians here — if you prefer degrees we can convert)
    lines.append("[ angles ]")
    lines.append(";  ai    aj    ak  funct     th0(deg)        k(kJ/mol/rad^2)")
    for (i, j, k_, th0_rad, k_rad, sigma_rad) in angles_fit:
        th0_deg = th0_rad * RAD2DEG
        lines.append(f"{i+1:5d} {j+1:5d} {k_+1:5d} {1:5d} {th0_deg:12.5f} {k_rad:16.2f}")
    lines.append("")

    lines.append("[ dihedrals ]")
    lines.append(";  ai    aj    ak    al  funct   phi0(deg)        k(kJ/mol)   mult")
    for (i,j,k,l,phi0_rad,kk,n) in dihedrals_fit:
        phi0_deg = phi0_rad * RAD2DEG
        lines.append(f"{i+1:5d} {j+1:5d} {k+1:5d} {l+1:5d} {1:5d} {phi0_deg:12.5f} {kk:12.3f} {n:6d}")

    Path(out_itp).write_text("\n".join(lines) + "\n")

# -----------------------------
# 5) The one function you call
# -----------------------------

def build_cg_from_xtb(
    smiles: str,
    final: List[str],
    mapping,
    ref_gro: str,
    ref_xtc: str,
    compound: str,
    out_prefix: str,
    T: float = 300.0,
):
    """
    End-to-end:
      final+mapping -> bead groups + bead types + bead bonds (reuses your code)
      smiles -> add attached H indices (RDKit)
      ref.gro/ref.xtc -> cg_ref.gro/cg_ref.xtc (center of geometry, with H)
      fit bonds+angles from cg_ref.xtc
      write report + itp

    Outputs:
      {out_prefix}_cg_ref.gro
      {out_prefix}_cg_ref.xtc
      {out_prefix}_bonded_report.txt
      {out_prefix}.itp
    """
    # bead grouping + types (heavy atoms only)
    bead_heavy_atoms, bead_types_raw = group_beads_by_type(final)  # :contentReference[oaicite:2]{index=2}
    bead_types_fixed = fix_beadtypes(bead_types_raw)               # :contentReference[oaicite:3]{index=3}

    # bead connectivity from mapping
    origin, connected, bond_types = compute_beads_connections(bead_heavy_atoms, mapping)  # :contentReference[oaicite:4]{index=4}
    n_beads = len(bead_heavy_atoms)

    # unique bonds list (0-based bead indices)
    bonds = sorted({(min(i, j), max(i, j)) for i, j in zip(origin, connected)})

    # map AA -> CG trajectory
    out_cg_gro = f"{out_prefix}_cg_ref.gro"
    out_cg_xtc = f"{out_prefix}_cg_ref.xtc"
    bead_names = [f"C{i+1}" for i in range(n_beads)]
    # Load AA once (also lets us fail early if ref.gro is corrupt)
    aa = md.load(ref_xtc, top=ref_gro)
    
    # Strict: assign H by true bond connectivity from SMILES (RDKit)
    heavy_to_h = build_heavy_to_hydrogens(smiles)
    
    # Expand each bead's atom list to include its attached H
    bead_atom_indices = expand_beads_with_hydrogens(bead_heavy_atoms, heavy_to_h)
    
    # Now map AA -> CG using the already-loaded aa
    cg = map_aa_to_cg_traj_from_loaded_aa(
        aa=aa,
        bead_atom_indices=bead_atom_indices,
        bead_names=bead_names,
        out_cg_gro=out_cg_gro,
        out_cg_xtc=out_cg_xtc,
    )

    # generate angles from bond graph + fit
    angles = generate_angle_list_from_bonds(n_beads, bonds)
    bonds_fit, angles_fit = fit_bonds_and_angles(cg, bonds, angles, T=T)
    
    dihedrals = generate_dihedral_list_from_bonds(n_beads, bonds)
    dihedrals_fit = fit_dihedrals_simple(cg, dihedrals, T=T)

    out_itp = f"{out_prefix}.itp"
    write_itp_with_bonds_angles_dihedrals(compound, bead_types_fixed, bonds_fit, angles_fit, dihedrals_fit, out_itp)

    return {
        "cg_gro": out_cg_gro,
        "cg_xtc": out_cg_xtc,
        "itp": out_itp,
        "n_beads": n_beads,
        "n_frames": cg.n_frames,
        "bonds": bonds,
        "angles": angles,
        "dihedrals": dihedrals
    }
    
